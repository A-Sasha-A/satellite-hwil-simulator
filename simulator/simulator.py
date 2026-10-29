# Import necessary libraries
import numpy as np
import itertools
import json
import time
import signal
from datetime import datetime, timedelta
from pathlib import Path
import serial
import sys
import multiprocessing as mp

# --- ASCII dashboard support (works both in Jupyter and terminal) ---
try:
    from IPython.display import clear_output
    _HAS_IPYTHON = True
except ImportError:
    _HAS_IPYTHON = False
    def clear_output(wait=True):
        pass

# --- Cartopy support (for the map, handled in the child process) ---
try:
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    _HAS_CARTOPY = True
except ImportError:
    _HAS_CARTOPY = False

# --- Map process module (must be in map_process.py in the same folder) ---
try:
    from map_process import map_process as _map_process_fn
    _HAS_MAP_PROC = True
except ImportError:
    _HAS_MAP_PROC = False


# INITIAL ORBIT CHARACTERISTICS (CUSTOMIZABLE)
semi_major_axis = 6878.14    # [km]
eccentricity = 0.0           # circular orbit (Sentinel-2 like): drag keeps it ~0;
                             # the firmware will only intervene if it grows due to anomalies
inclination = 97.4           # [deg]
raan = 178.0                 # [deg] LTAN approx 10:30 (daytime SSO, Sentinel-2-like)
argument_of_perigee = 0.0    # [deg]
true_anomaly = 0.0           # [deg]


# INITIAL TIME (CUSTOMIZABLE)
year = 2026
month = 10
day = 15
hour = 10
minute = 0
second = 0


# SATELLITE CONFIGURATION (CUSTOMIZABLE)
# NOTE: the 6 cartesian variables are only used if USE_KEPLERIAN_ELEMENTS = False
initial_pos_x = 6878.14   # [km]
initial_pos_y = 0.0       # [km]
initial_pos_z = 0.0       # [km]
initial_vel_x = 0.0       # [km/s]
initial_vel_y = 7.61      # [km/s]
initial_vel_z = 0.0       # [km/s]

# Attitude: body angles (3-2-1) relative to the local orbital frame LVLH
# (x = velocity, y = -orbit normal, z = nadir = camera face).
initial_roll = 0.0    # [deg]
initial_pitch = 0.0   # [deg]
initial_yaw = 0.0     # [deg]

# Initial body angular velocity, RELATIVE to the LVLH frame
# (if INITIAL_OMEGA_NADIR = True it is added to the orbital angular velocity required
#  to maintain nadir pointing; if False it represents the absolute inertial angular velocity).
initial_omega_x = 0.0   # [rad/s]
initial_omega_y = 0.0   # [rad/s]
initial_omega_z = 0.0   # [rad/s]

# Initial angular momentum of the 4 reaction wheels (along their own axis)
initial_wheel_momentum = np.zeros(4)   # [N m s]


# SATELLITE CHARACTERISTICS (CUSTOMIZABLE)
# Body axes: x = length, y = width, z = height. Solar panels extend along ±y.
# Camera face: +z (nadir). Rear face (Hall thruster): -x.
length = 0.8    # [m]
width = 0.8     # [m]
height = 1.2    # [m]
solar_panel_length = 1.5   # [m] (y-axis extension, per panel)
solar_panel_width = 0.6    # [m] (x-axis extension)
dry_mass = 120             # [kg]
xenon_mass = 15            # [kg]
chemical_propellant_mass = 15  # [kg]
solar_panel_mass = 2.0     # [kg] mass of each panel (ASSUMPTION, needed for the inertia tensor)
center_of_mass_offset = np.array([0.0, 0.0, 0.0])   # [m] CoM offset from the geometric center (body axes)


# =====================================================================
# ACTUATORS (CUSTOMIZABLE) - typical values for a ~150 kg microsatellite
# =====================================================================
g0 = 9.80665   # [m/s^2]

# --- Hall effect thruster (ONLY one, variable intensity) ---
# Mounted at the center of the rear face (-x), thrust along body +x = velocity direction (LVLH).
max_hall_thrust = 0.015            # [N]   (15 mN, 300-400 W class)
min_hall_fraction = 0.2            # [-]   below this command fraction the thruster shuts down (stability limit)
isp_hall = 1400.0                  # [s]
hall_position = np.array([-length / 2.0, 0.0, 0.0])   # [m] from geometric center
hall_thrust_direction = np.array([1.0, 0.0, 0.0])     # [-] force unit vector on the satellite (body frame)

# --- 8 chemical RCS thrusters (on/off boolean), fixed, one per corner ---
# Each thruster is at corner (sx, sy, sz) and expels gas in direction (sx, 0, sz)/sqrt(2),
# i.e., outwards and angled at 45° in the x-z plane (never towards the ±y solar panels).
# This layout provides ALL pure torques (±x, ±y, ±z) and pure translations ±x, ±z.
rcs_thrust = 0.5                     # [N] per thruster
isp_rcs = 210.0                      # [s] (monopropellant)

# --- 4 reaction wheels in a pyramidal configuration ---
# Axes: a_i = [cos(beta) cos(phi_i), cos(beta) sin(phi_i), sin(beta)], phi = 0, 90, 180, 270 deg.
# beta = 35.26° makes the torque/momentum capacity ISOTROPIC along the 3 body axes.
wheel_beta_angle = 35.26             # [deg] axis elevation relative to the x-y plane
max_wheel_torque = 0.020             # [N m] (20 mNm)
max_wheel_momentum = 0.5             # [N m s] (0.5 Nms)
max_wheel_speed = 6000.0 * 2.0 * np.pi / 60.0   # [rad/s] (6000 rpm)
wheel_inertia = max_wheel_momentum / max_wheel_speed   # [kg m^2] consistent with h_max = J * Omega_max

# --- 3 orthogonal magnetorquers (variable intensity), one per body axis ---
max_mtq_dipole = 20.0                # [A m^2]

# Wheels are commanded in SPEED [RPM, with sign] by the ESP32
WHEEL_MODE = "SPEED"                 # "SPEED": speed setpoint (RPM); "TORQUE": normalized torque as before
wheel_loop_tau = 2.0                 # [s] speed loop time constant of the wheel driver
                                     #     (keep dt_sim <= 1 s; must match WHEEL_LOOP_TAU in the firmware)
wheel_loop_gain = wheel_inertia / wheel_loop_tau   # [N m s] torque = gain * (Omega_cmd - Omega)


# MODELING OPTIONS (CUSTOMIZABLE)
USE_KEPLERIAN_ELEMENTS = True      # True: initial state from orbital elements; False: from cartesian values
# "FREE": attitude propagated with Euler's equations (wheels, magnetorquers, RCS affect the attitude)
# "LVLH": ideal prescribed nadir-pointing attitude (attitude actuators have NO effect; thrust does)
ATTITUDE_MODE = "FREE"
INITIAL_OMEGA_NADIR = True         # True: starts with orbital angular velocity (nadir pointing)
PANELS_FOLLOW_SUN = True           # True: panels rotate around y to point at the Sun
Cd = 2.2                           # [-] aerodynamic drag coefficient
solar_activity_factor = 1.0        # [-] density multiplier (1 = mean activity; >1 high activity)
# Optical coefficients (ASSUMPTIONS): rho_s = specular reflection, rho_d = diffuse reflection
RHO_S_BUS, RHO_D_BUS = 0.10, 0.20
RHO_S_PANELS, RHO_D_PANELS = 0.05, 0.10
# Enabling individual contributions (useful for validating the model)
PERTURBATIONS = {"J2": True, "SUN": True, "MOON": True, "DRAG": True, "SRP": True}

# =====================================================================
# SPACE DEBRIS — model parameters
# =====================================================================
DEBRIS_ACTIVE            = True
DEBRIS_ACCELERATION      = 1.0           # [-] 1.0 = real physics (impacts are extremely rare, see NOTES)
                                         #     1e6 = mild demo (~3 impacts per simulated orbit)
                                         #     1e7 = visible demo (~30 impacts)
                                         #     1e8 = aggressive demo (~300 impacts)
DAMAGE_PER_IMPACT        = 1.0e-3        # [-] fraction of integrity lost per impact (0.1%)
INITIAL_PANEL_INTEGRITY  = 1.0           # [-] 1.0 = 100% Beginning Of Life
SOLAR_PANEL_POWER_BOL_W  = 270.0         # [W] electrical power at BOL (proxy: 1.8 m² × 150 W/m²)
DEBRIS_SEED              = 42            # [-] None = random; int = reproducible across runs

# --- Large debris (>10 cm) and collision avoidance maneuvers ---
LARGE_DEBRIS_ACTIVE = True
LARGE_DEBRIS_ACCELERATION = 1.0e5        # 1.0 = real physics (rare events); 1e10 = ~1 event every 20 sim minutes
LARGE_DEBRIS_REL_V_KM_S = 10.0           # typical relative velocity [km/s]
CAM_WARNING_TIME_S = 3600.0              # collision avoidance warning time (1 hour)
CAM_ALTITUDE_DELTA_KM = 2.0              # required altitude increase to avoid the debris

# Large debris density table [debris / km³]
_LARGE_DEBRIS_ALTITUDES = np.array([200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 2000], dtype=float)
_LARGE_DEBRIS_DENSITY = np.array([
    1.0e-10, 2.0e-10, 4.0e-10, 8.0e-10, 1.5e-9,
    2.0e-9, 3.0e-9, 2.5e-9, 1.5e-9, 1.0e-9,
    1.2e-9, 1.8e-9, 2.0e-9, 1.0e-9, 5.0e-10, 2.0e-10, 5.0e-11
])

def get_large_debris_density(altitude_km):
    """Density [debris/km³] for objects >10 cm. Returns 0 outside the table range."""
    if altitude_km < _LARGE_DEBRIS_ALTITUDES[0] or altitude_km > _LARGE_DEBRIS_ALTITUDES[-1]:
        return 0.0
    return float(np.interp(altitude_km, _LARGE_DEBRIS_ALTITUDES, _LARGE_DEBRIS_DENSITY))


# CELESTIAL BODIES CHARACTERISTICS
sun_mass = 1.989e30                  # [kg]
mu_sun = 1.32712440018e11            # [km^3/s^2]
sun_earth_distance = 149600000.0     # [km]
sun_radius = 696000.0                # [km]
sun_radiation_emission = 1361.0      # [W/m^2]
moon_mass = 7.342e22                 # [kg]
mu_moon = 4902.800066                # [km^3/s^2]
earth_moon_distance = 384400.0       # [km]
moon_reflectivity = 0.12             # [-] (unused: would be needed for moon-reflected radiation pressure / thermal)
solar_radiation_pressure = 4.56e-6   # [N/m^2] (informational: the code uses sun_radiation_emission / c)
earth_albedo = 0.30                  # [-] (unused: would be needed for Earth albedo / IR)


# SUPPORT VARIABLES
universal_gravitational_constant = 6.6743e-11   # [N*(m^2)/(kg^2)] (unused, kept for reference)
earth_mass = 5.97219e24                         # [kg]
earth_equatorial_radius = 6378.14               # [km]
J2 = 0.00108263                                 # [-] J2 zonal coefficient (previously called Earth flattening coefficient)
earth_flattening = 1.0 / 298.257223563          # [-] true geometric flattening (WGS84)
earth_rotation_rate = 7.292115e-5               # [rad/s]
mu = 398600.4418                                # [km^3/s^2]
AU_KM = 149597870.7                             # [km] astronomical unit
speed_of_light = 299792.458                     # [km/s]
# Geomagnetic field: tilted dipole from IGRF-13 (2020.0 epoch, coefficients in nT, secular variation in nT/year)
igrf_reference_radius = 6371.2                  # [km]
IGRF_G10, IGRF_G11, IGRF_H11 = -29404.8, -1450.9, 4652.5
IGRF_SV_G10, IGRF_SV_G11, IGRF_SV_H11 = 5.7, 7.4, -25.9


# =====================================================================
# SPACE DEBRIS — solar panel impact flux
# Based on REAL FLUX [impacts / m² / year] for debris > 1 mm.
# Anchors: ESA Space Debris Mitigation Handbook (Tab. 1.3-1) and MASTER-8.
# Intermediate values are interpolated between known altitudes.
# =====================================================================
_FLUX_ALTITUDES_KM = np.array([200, 400, 600, 800, 1000, 1200, 1500, 2000], dtype=float)

# Flux [impacts / m² / year] for debris > 1 mm (diameter 1 mm – 1 cm)
_IMPACT_FLUX_M2_YEAR = np.array([
    1.0e-3,   # 200 km  (atmosphere sweeps debris away)
    2.6e-3,   # 400 km  (ESA: 1 impact every 3.9 years on 100 m²)
    5.0e-3,   # 600 km  (interpolated)
    1.0e-2,   # 800 km  (ESA: 1 impact every 1.0 year on 100 m²)
    8.0e-3,   # 1000 km (just below the peak)
    5.0e-3,   # 1200 km
    2.0e-3,   # 1500 km
    5.0e-4,   # 2000 km
])


def get_impact_flux(altitude_km):
    """Flux [impacts / m² / year] for debris >1 mm. Returns 0 outside table range."""
    if altitude_km < _FLUX_ALTITUDES_KM[0] or altitude_km > _FLUX_ALTITUDES_KM[-1]:
        return 0.0
    return float(np.interp(altitude_km, _FLUX_ALTITUDES_KM, _IMPACT_FLUX_M2_YEAR))


def update_panel_integrity(integrity, total_impacts, altitude_km,
                           panel_area_m2, dt_seconds, rng):
    """
    One step of panel degradation due to > 1 mm debris impacts.

    Physics:
      λ (expected impacts in dt) = flux[impacts/m²/year] × area[m²] × dt[s] / seconds_in_a_year
      impacts_dt ~ Poisson(λ × DEBRIS_ACCELERATION)
      integrity *= (1 - DAMAGE_PER_IMPACT)^impacts_dt

    Returns (integrity, total_impacts, impacts_dt, real_lambda).
    `real_lambda` is the physical rate (NOT accelerated), useful for the dashboard.
    """
    if not DEBRIS_ACTIVE:
        return integrity, total_impacts, 0, 0.0

    annual_flux = get_impact_flux(altitude_km)
    if annual_flux <= 0.0 or panel_area_m2 <= 0.0:
        return integrity, total_impacts, 0, 0.0

    seconds_in_a_year = 365.25 * 24.0 * 3600.0
    real_lambda = annual_flux * panel_area_m2 * (dt_seconds / seconds_in_a_year)
    lam = real_lambda * DEBRIS_ACCELERATION

    # Poisson draw. If lam is tiny, this will almost always be 0, which is exactly what we want.
    impacts_dt = int(rng.poisson(lam))

    if impacts_dt > 0:
        integrity = integrity * (1.0 - DAMAGE_PER_IMPACT) ** impacts_dt
        integrity = max(0.0, integrity)
        total_impacts += impacts_dt

    return integrity, total_impacts, impacts_dt, real_lambda


# =====================================================================
# UTILITIES: KEPLERIAN ELEMENTS -> ECI STATE
# =====================================================================
def keplerian_to_eci(a, e, i_deg, raan_deg, argp_deg, nu_deg, mu):
    """Converts classical orbital elements to ECI position [km] and velocity [km/s]."""
    # Standard Vallado formulas. Nothing fancy here, just trig.
    i, O, w, nu = np.radians([i_deg, raan_deg, argp_deg, nu_deg])
    p = a * (1 - e ** 2)
    r_mag = p / (1 + e * np.cos(nu))
    
    # Position and velocity in the perifocal frame (PQW)
    r_pf = r_mag * np.array([np.cos(nu), np.sin(nu), 0.0])
    v_pf = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

    # Rotation matrix from PQW to ECI
    cO, sO, cw, sw, ci, si = np.cos(O), np.sin(O), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    R = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci,  sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si,                 cw * si,                  ci]
    ])
    
    return R @ r_pf, R @ v_pf

# =====================================================================
# UTILITIES: ROTATIONS, QUATERNIONS, ATTITUDE, INERTIA
# =====================================================================
def rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def rot_y(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def quat_multiply(p, q):
    """Quaternion product [w, x, y, z]. 
    Quaternion math is always a headache, but at least we don't have to worry about gimbal lock."""
    pw, px, py, pz = p
    qw, qx, qy, qz = q
    return np.array([
        pw * qw - px * qx - py * qy - pz * qz,
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw
    ])


def quat_to_matrix(q):
    """Rotation matrix body -> ECI from quaternion [w, x, y, z]."""
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w),     2 * (x * z + y * w)],
        [2 * (x * y + z * w),     1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w),     2 * (y * z + x * w),     1 - 2 * (x * x + y * y)]
    ])


def matrix_to_quat(R):
    """Quaternion [w, x, y, z] from rotation matrix (Shepperd's method). 
    Numerically stable, just trust the process on the branches."""
    tr = np.trace(R)
    if tr > 0:
        s = np.sqrt(tr + 1.0) * 2.0
        w = 0.25 * s
        x = (R[2, 1] - R[1, 2]) / s
        y = (R[0, 2] - R[2, 0]) / s
        z = (R[1, 0] - R[0, 1]) / s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2.0
        w = (R[2, 1] - R[1, 2]) / s
        x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s
        z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2.0
        w = (R[0, 2] - R[2, 0]) / s
        x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s
        z = (R[1, 2] + R[2, 1]) / s
    else:
        s = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2.0
        w = (R[1, 0] - R[0, 1]) / s
        x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s
        z = 0.25 * s
    q = np.array([w, x, y, z])
    return q / np.linalg.norm(q)


def lvlh_matrix(r_vec, v_vec):
    """LVLH -> ECI matrix. Columns: x ~ velocity, y = -orbit normal, z = nadir."""
    # Nadir pointing is just the negative position vector
    z = -r_vec / np.linalg.norm(r_vec)
    h = np.cross(r_vec, v_vec) # orbit normal
    y = -h / np.linalg.norm(h)
    x = np.cross(y, z)         # completes the right-handed triad
    return np.column_stack([x, y, z])


# Constant body rotation relative to LVLH (3-2-1 sequence: yaw, pitch, roll)
ATTITUDE_OFFSET_MATRIX = (rot_z(np.radians(initial_yaw))
                          @ rot_y(np.radians(initial_pitch))
                          @ rot_x(np.radians(initial_roll)))


def attitude_matrix(r_vec, v_vec, q):
    """Body -> ECI matrix according to the selected attitude mode."""
    if ATTITUDE_MODE == "LVLH":
        return lvlh_matrix(r_vec, v_vec) @ ATTITUDE_OFFSET_MATRIX
    return quat_to_matrix(q)


def calculate_inertia_tensor(m_total=None):
    """
    Inertia tensor [kg m^2] in the body frame, for the instantaneous total mass m_total [kg]:
    bus = homogeneous rectangular box (residual propellant is assumed distributed inside the bus),
    panels = thin plates at ±y (neglecting booms/actuator masses).
    """
    if m_total is None:
        m_total = dry_mass + xenon_mass + chemical_propellant_mass
        
    # We remove the panels' mass from the main bus mass
    m_bus = m_total - 2.0 * solar_panel_mass
    
    # Box inertia for the main bus
    Ixx = m_bus * (width ** 2 + height ** 2) / 12.0
    Iyy = m_bus * (length ** 2 + height ** 2) / 12.0
    Izz = m_bus * (length ** 2 + width ** 2) / 12.0
    
    # Distance from bus center to panel center along y
    d = width / 2.0 + solar_panel_length / 2.0
    
    # Add the panels' contribution (parallel axis theorem)
    mp, Lp, Wp = solar_panel_mass, solar_panel_length, solar_panel_width
    Ixx += 2.0 * (mp * Lp ** 2 / 12.0 + mp * d ** 2)
    Iyy += 2.0 * (mp * Wp ** 2 / 12.0)
    Izz += 2.0 * (mp * (Wp ** 2 + Lp ** 2) / 12.0 + mp * d ** 2)
    
    # The satellite is roughly symmetric, so the products of inertia are zero
    return np.diag([Ixx, Iyy, Izz])


# =====================================================================
# ACTUATOR GEOMETRY (constant matrices in the body frame)
# =====================================================================
# Reaction wheels: matrix A (3 x 4) with the axes of the 4 wheels as columns
_beta_rad = np.radians(wheel_beta_angle)
_phi_rad = np.radians([0.0, 90.0, 180.0, 270.0])
WHEEL_AXES = np.array([np.cos(_beta_rad) * np.cos(_phi_rad),
                       np.cos(_beta_rad) * np.sin(_phi_rad),
                       np.sin(_beta_rad) * np.ones(4)])
WHEEL_AXES_PINV = np.linalg.pinv(WHEEL_AXES)   # 4 x 3, minimum norm distribution

# RCS: 8 thrusters, one per corner, force on the satellite = opposite to the jet
_CORNERS = [(sx, sy, sz) for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]
RCS_POSITIONS = np.array([[sx * length / 2.0, sy * width / 2.0, sz * height / 2.0]
                          for sx, sy, sz in _CORNERS]).T                       # 3 x 8 [m]
RCS_DIRECTIONS = np.array([-np.array([sx, 0.0, sz]) / np.sqrt(2.0)
                           for sx, sy, sz in _CORNERS]).T                      # 3 x 8 (force unit vector)
# Force [N] and torque [N m] of each active thruster (torque relative to CoM)
UNIT_RCS_FORCE = rcs_thrust * RCS_DIRECTIONS
UNIT_RCS_TORQUE = np.array([np.cross(RCS_POSITIONS[:, i] - center_of_mass_offset, UNIT_RCS_FORCE[:, i])
                            for i in range(8)]).T

HALL_DIRECTION = hall_thrust_direction / np.linalg.norm(hall_thrust_direction)


def calculate_rcs_combinations():
    """
    Enumerates the 256 on/off RCS combinations and returns, for each pure motion,
    the combination (array of 8 booleans) with maximum effect. Useful as a lookup table for the ESP32.
    Pure motions: torque around ±x, ±y, ±z with zero net force; force along ±x, ±z with zero net torque.
    """
    # Brute force over all 256 possible states. It's fast enough at startup.
    combinations = np.array(list(itertools.product([0, 1], repeat=8)), dtype=float)
    F = combinations @ UNIT_RCS_FORCE.T
    M = combinations @ UNIT_RCS_TORQUE.T
    tol = 1e-9
    results = {}
    axes = ["x", "y", "z"]
    
    # Pure torques
    for ax in range(3):
        others = [j for j in range(3) if j != ax]
        for s, sign in ((1, "+"), (-1, "-")):
            ok = (np.abs(F).max(1) < tol) & (np.abs(M[:, others]).max(1) < tol) & (s * M[:, ax] > tol)
            if ok.any():
                idx = np.where(ok)[0]
                k = idx[np.argmax(s * M[idx, ax])]
                results[f"torque_{sign}{axes[ax]}"] = combinations[k].astype(bool)
                
    # Pure forces (only x and z are fully decoupled by this geometry)
    for ax in (0, 2):
        others = [j for j in range(3) if j != ax]
        for s, sign in ((1, "+"), (-1, "-")):
            ok = (np.abs(M).max(1) < tol) & (np.abs(F[:, others]).max(1) < tol) & (s * F[:, ax] > tol)
            if ok.any():
                idx = np.where(ok)[0]
                k = idx[np.argmax(s * F[idx, ax])]
                results[f"force_{sign}{axes[ax]}"] = combinations[k].astype(bool)
                
    return results


RCS_COMBINATIONS = calculate_rcs_combinations()


# =====================================================================
# SATELLITE GEOMETRY FROM ATTITUDE: normals, areas, optical coefficients
# =====================================================================
def satellite_surfaces(R_b2i, sun_dir_eci, length, width, height,
                       panel_length, panel_width):
    """
    Returns (normals_eci [N x 3], areas [m^2], rho_s, rho_d) of all flat surfaces:
    6 bus faces + 2 panel sides (the two panels are assumed parallel).
    Normals are outgoing. Mutual shading between panels and bus is not considered.
    """
    A_x = width * height          # faces with normal ±x
    A_y = length * height         # faces with normal ±y
    A_z = length * width          # faces with normal ±z
    A_panels = 2.0 * panel_length * panel_width   # area of both panels (one side)

    # Panel normal in the body frame (they rotate around y)
    if PANELS_FOLLOW_SUN:
        s_b = R_b2i.T @ sun_dir_eci
        n_p = np.array([s_b[0], 0.0, s_b[2]])
        norm = np.linalg.norm(n_p)
        n_p = n_p / norm if norm > 1e-9 else np.array([0.0, 0.0, 1.0])
    else:
        n_p = np.array([0.0, 0.0, 1.0])

    normals_body = np.array([
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1],
        n_p, -n_p
    ], dtype=float)
    
    areas = np.array([A_x, A_x, A_y, A_y, A_z, A_z, A_panels, A_panels])
    rho_s = np.array([RHO_S_BUS] * 6 + [RHO_S_PANELS] * 2)
    rho_d = np.array([RHO_D_BUS] * 6 + [RHO_D_PANELS] * 2)
    normals_eci = (R_b2i @ normals_body.T).T
    return normals_eci, areas, rho_s, rho_d


# =====================================================================
# ATMOSPHERE: GEODETIC ALTITUDE AND ALTITUDE-DEPENDENT DENSITY
# =====================================================================
def geodetic_altitude(r_vec):
    """Altitude [km] on the WGS84 ellipsoid (iterative); takes into account Earth's flattening."""
    # Standard iterative algorithm. Converges in a few iterations.
    a = earth_equatorial_radius
    e2 = earth_flattening * (2.0 - earth_flattening)
    x, y, z = r_vec
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1.0 - e2))
    h = 0.0
    for _ in range(6):
        N = a / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
        if abs(lat) < 1.3:
            h = p / np.cos(lat) - N
        else:
            h = z / np.sin(lat) - N * (1.0 - e2)
        lat = np.arctan2(z, p * (1.0 - e2 * N / (N + h)))
    return h


# Exponential layered model (Vallado, reference atmosphere for mean solar activity)
# columns: base altitude h0 [km], base density rho0 [kg/m^3], scale height H [km]
_ATMOSPHERE_TABLE = np.array([
    [0,    1.225,     7.249],   [25,   3.899e-2, 6.349],   [30,   1.774e-2, 6.682],
    [40,   3.972e-3,  7.554],   [50,   1.057e-3, 8.382],   [60,   3.206e-4, 7.714],
    [70,   8.770e-5,  6.549],   [80,   1.905e-5, 5.799],   [90,   3.396e-6, 5.382],
    [100,  5.297e-7,  5.877],   [110,  9.661e-8, 7.263],   [120,  2.438e-8, 9.473],
    [130,  8.484e-9,  12.636],  [140,  3.845e-9, 16.149],  [150,  2.070e-9, 22.523],
    [180,  5.464e-10, 29.740],  [200,  2.789e-10, 37.105], [250,  7.248e-11, 45.546],
    [300,  2.418e-11, 53.628],  [350,  9.518e-12, 53.298], [400,  3.725e-12, 58.515],
    [450,  1.585e-12, 60.828],  [500,  6.967e-13, 63.822], [600,  1.454e-13, 71.835],
    [700,  3.614e-14, 88.667],  [800,  1.170e-14, 124.64], [900,  5.245e-15, 181.05],
    [1000, 3.019e-15, 268.00]
])
_H0_TABLE = _ATMOSPHERE_TABLE[:, 0]


def atmospheric_density(r_vec):
    """Atmospheric density [kg/m^3] as a function of geodetic altitude (continuous, no hard cuts)."""
    h = max(geodetic_altitude(r_vec), 0.0)
    
    # Find the right layer in the table
    idx = int(np.searchsorted(_H0_TABLE, h, side="right") - 1)
    idx = min(max(idx, 0), len(_H0_TABLE) - 1)
    h0, rho0, H = _ATMOSPHERE_TABLE[idx]
    
    return solar_activity_factor * rho0 * np.exp(-(h - h0) / H)


# =====================================================================
# GEOMAGNETIC FIELD (tilted IGRF-13 dipole, needed by the magnetorquers)
# =====================================================================
def geomagnetic_field_eci(r_vec, jd):
    """
    Earth's magnetic field B [T] in the ECI frame at position r_vec [km] and Julian day jd.
    Tilted dipole model (IGRF-13 degree 1 with secular variation): typical accuracy 10-30%
    in LEO. For higher fidelity, this could be replaced by a full IGRF (degree 13) model.
    """
    # Calculate how many years have passed since the 2020.0 epoch
    years = (jd - 2451545.0) / 365.25 + 2000.0 - 2020.0
    g10 = IGRF_G10 + IGRF_SV_G10 * years
    g11 = IGRF_G11 + IGRF_SV_G11 * years
    h11 = IGRF_H11 + IGRF_SV_H11 * years

    # Earth rotation angle (GMST) and dipole moment ECEF -> ECI
    gmst = np.radians((280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360.0)
    M_ecef = np.array([g11, h11, g10])      # [nT]
    M_eci = rot_z(gmst) @ M_ecef

    r = np.linalg.norm(r_vec)
    r_hat = r_vec / r
    
    # Standard tilted dipole equation
    B_nT = (igrf_reference_radius / r) ** 3 * (3.0 * np.dot(M_eci, r_hat) * r_hat - M_eci)
    
    # Convert from nT to Tesla
    return B_nT * 1e-9


# =====================================================================
# CELESTIAL BODY POSITIONS (Earth, Sun, Moon, Satellite)
# Rigorous astrodynamics version using Julian Time and rotation matrices.
# It's a bit of a beast, but it does the job.
# =====================================================================
def calculate_body_positions(dt_obj, r_satellite_eci):
    """
    Returns the positions in the Earth-Centered Inertial (ECI) frame [km] of:
    - Earth (at the origin: [0, 0, 0])
    - Sun (Analytical model based on Julian Time and J2000.0)
    - Moon (Elliptical model with equation of center, ECI rotation)
    - Satellite (passed as an ECI parameter)
    Time is passed as a datetime object (microseconds included).
    """
    year_, month_, day_ = dt_obj.year, dt_obj.month, dt_obj.day
    seconds_ = dt_obj.second + dt_obj.microsecond * 1e-6

    # 1. Calculate Julian Date (JD) and Julian Centuries (T) from J2000.0
    if month_ <= 2:
        a_mod = year_ - 1
        m_mod = month_ + 12
    else:
        a_mod = year_
        m_mod = month_

    A = int(a_mod / 100)
    B = 2 - A + int(A / 4)
    jd = int(365.25 * (a_mod + 4716)) + int(30.6001 * (m_mod + 1)) + day_ + B - 1524.5
    day_fraction = (dt_obj.hour + dt_obj.minute / 60.0 + seconds_ / 3600.0) / 24.0
    jd_total = jd + day_fraction

    # Julian centuries since J2000.0 (Jan 1st 2000, 12:00 TT)
    T = (jd_total - 2451545.0) / 36525.0
    d = jd_total - 2451545.0  # Days since J2000.0

    # --- 2. SUN POSITION (High-precision analytical model) ---
    # Standard low-precision solar ephemeris. Good enough for our LEO drag/SRP calculations.
    L0 = 280.46646 + 36000.76983 * T
    L0 = L0 % 360.0

    M_sun = 357.52911 + 35999.05029 * T
    M_sun_rad = np.radians(M_sun % 360.0)

    C_sun = (1.914602 - 0.004817 * T - 0.000014 * T ** 2) * np.sin(M_sun_rad) + \
            (0.019993 - 0.000101 * T) * np.sin(2.0 * M_sun_rad) + \
            0.000289 * np.sin(3.0 * M_sun_rad)

    sun_true_long = L0 + C_sun
    sun_true_long_rad = np.radians(sun_true_long % 360.0)

    e_earth = 0.016708634 - 0.000042037 * T - 0.0000001267 * T ** 2
    true_anomaly_sun = M_sun + C_sun
    true_anomaly_sun_rad = np.radians(true_anomaly_sun % 360.0)

    dist_au = (1.000001018 * (1.0 - e_earth ** 2)) / (1.0 + e_earth * np.cos(true_anomaly_sun_rad))
    dist_km = dist_au * AU_KM

    obliquity = 23.439291 - 0.0130042 * T - 0.00000016 * T ** 2 + 0.000000504 * T ** 3
    eps_rad = np.radians(obliquity)

    x_ecl = dist_km * np.cos(sun_true_long_rad)
    y_ecl = dist_km * np.sin(sun_true_long_rad)
    z_ecl = 0.0

    # Rotate from ecliptic to ECI equatorial frame
    x_sun_eci = x_ecl
    y_sun_eci = y_ecl * np.cos(eps_rad) - z_ecl * np.sin(eps_rad)
    z_sun_eci = y_ecl * np.sin(eps_rad) + z_ecl * np.cos(eps_rad)
    pos_sun = np.array([x_sun_eci, y_sun_eci, z_sun_eci])  # [km]

    # --- 3. MOON POSITION (Elliptical model, equation of center) ---
    # Don't ask me to derive this from scratch, it's a simplified version of the standard Meeus formulas.
    # We skip the Evection and Variation terms, they are too small for our purposes.
    L_moon = 218.316447 + 13.17639648 * d
    L_moon_rad = np.radians(L_moon % 360.0)

    M_moon = 134.963396 + 13.06499295 * d
    M_moon_rad = np.radians(M_moon % 360.0)

    F_moon = 93.272095 + 13.22935070 * d
    F_moon_rad = np.radians(F_moon % 360.0)

    # Longitude correction: just the equation of center (6.289° * sin M).
    l_pert = L_moon_rad + np.radians(6.289 * np.sin(M_moon_rad))
    lat_moon = np.radians(5.128 * np.sin(F_moon_rad))

    dist_moon = 385001.0 - 20905.0 * np.cos(M_moon_rad)

    x_moon_ecl = dist_moon * np.cos(lat_moon) * np.cos(l_pert)
    y_moon_ecl = dist_moon * np.cos(lat_moon) * np.sin(l_pert)
    z_moon_ecl = dist_moon * np.sin(lat_moon)

    # Ecliptic -> ECI
    x_moon_eci = x_moon_ecl
    y_moon_eci = y_moon_ecl * np.cos(eps_rad) - z_moon_ecl * np.sin(eps_rad)
    z_moon_eci = y_moon_ecl * np.sin(eps_rad) + z_moon_ecl * np.cos(eps_rad)
    pos_moon = np.array([x_moon_eci, y_moon_eci, z_moon_eci])  # [km]

    # --- 4. Earth and Satellite ---
    pos_earth = np.array([0.0, 0.0, 0.0])  # [km]
    pos_satellite = np.array(r_satellite_eci)  # [km]

    return {
        "Earth": pos_earth,
        "Sun": pos_sun,
        "Moon": pos_moon,
        "Satellite": pos_satellite,
        "JD": jd_total
    }


# =====================================================================
# ATMOSPHERIC DRAG ACCELERATION
# Variable density with altitude, effective area depends on attitude
# =====================================================================
def calculate_drag_acceleration(r_vec, v_vec, total_mass, length, width, height,
                                panel_length, panel_width, r_eq, omega_earth,
                                R_b2i, r_sun):
    """
    Calculates the vector acceleration due to atmospheric drag.
    - density: exponential layered model based on geodetic altitude
    - area: sum of the areas of the surfaces exposed to the flow (projection along v_rel)
    """
    rho = atmospheric_density(r_vec)

    # Atmosphere velocity due to Earth's rotation (rigid co-rotation) [km/s]
    v_atmosphere = np.cross(np.array([0.0, 0.0, omega_earth]), r_vec)

    v_rel = v_vec - v_atmosphere
    v_rel_mag = np.linalg.norm(v_rel)  # km/s

    if v_rel_mag == 0:
        return np.array([0.0, 0.0, 0.0])
    v_rel_unit = v_rel / v_rel_mag

    # Effective exposed area (m^2) based on attitude: sum up the surfaces facing the flow
    u_sun = r_sun - r_vec
    u_sun = u_sun / np.linalg.norm(u_sun)
    normals, areas, _, _ = satellite_surfaces(R_b2i, u_sun, length, width, height,
                                              panel_length, panel_width)
    cos_flow = normals @ v_rel_unit
    effective_area = float(np.sum(areas * np.clip(cos_flow, 0.0, None)))

    # Classic drag equation: F = 0.5 * rho * v^2 * Cd * A
    v_rel_ms = v_rel_mag * 1000.0
    acc_magnitude_ms2 = 0.5 * rho * (v_rel_ms ** 2) * Cd * (effective_area / total_mass)

    acc_magnitude_kms2 = acc_magnitude_ms2 / 1000.0
    a_drag = -acc_magnitude_kms2 * v_rel_unit

    return a_drag


# =====================================================================
# J2 PERTURBATION ACCELERATION
# Earth's oblateness (equatorial bulge) effect
# =====================================================================
def calculate_j2_acceleration(r_vec, mu, r_eq, J2):
    """
    Calculates the vector acceleration due to the J2 perturbation (Earth's flattening).
    r_vec: ECI position vector [x, y, z] in km
    mu: Standard gravitational parameter (km^3 / s^2)
    r_eq: Earth's equatorial radius (km)
    J2: Zonal coefficient J2 (dimensionless)
    """
    x, y, z = r_vec
    r = np.linalg.norm(r_vec)

    if r == 0.0:
        return np.array([0.0, 0.0, 0.0])

    r_sq = r ** 2
    z_sq = z ** 2

    # Standard J2 formula. Gets the biggest kick in when passing over the poles.
    factor = (1.5 * J2 * mu * (r_eq ** 2)) / (r ** 5)
    z_ratio_sq = z_sq / r_sq

    ax = factor * x * (5.0 * z_ratio_sq - 1.0)
    ay = factor * y * (5.0 * z_ratio_sq - 1.0)
    az = factor * z * (5.0 * z_ratio_sq - 3.0)

    return np.array([ax, ay, az])


# =====================================================================
# THIRD-BODY PERTURBATION ACCELERATION (Sun or Moon)
# =====================================================================
def calculate_third_body_acceleration(r_sat, r_body, mu_body):
    """
    Calculates the perturbing acceleration on the satellite due to a third body (Sun or Moon).
    r_sat: Satellite position vector relative to Earth (ECI) in km
    r_body: Celestial body position vector (Sun or Moon) relative to Earth (ECI) in km
    mu_body: Standard gravitational parameter of the third body (km^3 / s^2)
    """
    r_sat_body = r_sat - r_body

    d_sat_body = np.linalg.norm(r_sat_body)
    d_body_earth = np.linalg.norm(r_body)

    if d_sat_body == 0.0 or d_body_earth == 0.0:
        return np.array([0.0, 0.0, 0.0])

    # N-body perturbation formula (direct + indirect term)
    # Basically: the pull from the Sun on the satellite, minus the pull from the Sun on the Earth.
    term1 = r_sat_body / (d_sat_body ** 3)
    term2 = r_body / (d_body_earth ** 3)

    a_third_body = -mu_body * (term1 + term2)

    return a_third_body


# =====================================================================
# ILLUMINATION FACTOR: conical model with umbra and penumbra
# =====================================================================
def illumination_factor(r_sat, r_sun):
    """
    Fraction of the solar disk visible from the satellite (1 = full sun, 0 = total umbra).
    Conical model (Montenbruck & Gill): overlapping apparent disks of Sun and Earth.
    """
    r_sat_sun = r_sun - r_sat
    d_sun = np.linalg.norm(r_sat_sun)
    d_earth = np.linalg.norm(r_sat)

    a = np.arcsin(sun_radius / d_sun)             # Sun apparent angular radius
    b = np.arcsin(earth_equatorial_radius / d_earth)  # Earth apparent angular radius
    cos_c = np.dot(-r_sat, r_sat_sun) / (d_earth * d_sun)
    c = np.arccos(np.clip(cos_c, -1.0, 1.0))      # Angular separation of the centers

    if c >= a + b:
        return 1.0                                # No occultation
    if c <= b - a:
        return 0.0                                # Total umbra
    if c <= a - b:
        return 1.0 - (b ** 2) / (a ** 2)          # Annular eclipse (not for Earth, just for completeness)

    # Penumbra: partial overlap of the two disks
    x = (c ** 2 + a ** 2 - b ** 2) / (2.0 * c)
    y = np.sqrt(max(a ** 2 - x ** 2, 0.0))
    area = a ** 2 * np.arccos(np.clip(x / a, -1, 1)) + b ** 2 * np.arccos(np.clip((c - x) / b, -1, 1)) - c * y
    return 1.0 - area / (np.pi * a ** 2)


# =====================================================================
# SOLAR RADIATION PRESSURE (SRP) ACCELERATION
# Flat surface model with optical coefficients, penumbra and (1 AU / d)^2
# =====================================================================
def calculate_srp_acceleration(r_sat, r_sun, total_mass, length, width, height,
                               panel_length, panel_width, R_b2i):
    """
    Calculates the vector acceleration due to Solar Radiation Pressure (SRP).
    For each illuminated surface (normal n, cos(theta) = n . s):
      F = -P A cos(theta) [ (1 - rho_s) s + 2 (rho_s cos(theta) + rho_d / 3) n ]
    where s is the satellite -> Sun unit vector.
    """
    c_ms = speed_of_light * 1000.0   # m/s

    r_sat_sun = r_sun - r_sat
    d_sun_sat = np.linalg.norm(r_sat_sun)   # km
    if d_sun_sat == 0.0:
        return np.array([0.0, 0.0, 0.0])
    s_hat = r_sat_sun / d_sun_sat

    # Illumination (umbra + penumbra)
    illum = illumination_factor(r_sat, r_sun)
    if illum <= 0.0:
        return np.array([0.0, 0.0, 0.0])

    # Radiation pressure at the actual distance from the Sun [N/m^2]
    P = (sun_radiation_emission / c_ms) * (AU_KM / d_sun_sat) ** 2

    normals, areas, rho_s, rho_d = satellite_surfaces(R_b2i, s_hat, length, width, height,
                                                      panel_length, panel_width)
    cos_t = normals @ s_hat
    illuminated = cos_t > 0.0

    force = np.zeros(3)   # [N]
    for i in np.where(illuminated)[0]:
        ct = cos_t[i]
        # Classic SRP force equation. The specular and diffuse reflections matter a lot here.
        force += -P * areas[i] * ct * ((1.0 - rho_s[i]) * s_hat
                                       + 2.0 * (rho_s[i] * ct + rho_d[i] / 3.0) * normals[i])

    acc_ms2 = illum * force / total_mass
    return acc_ms2 / 1000.0   # [km/s^2]


# =====================================================================
# ACTUATORS: COMMANDS, FORCES AND TORQUES
# =====================================================================
def null_commands():
    """Idle commands: all actuators off."""
    return {"ion_thruster": 0.0, "wheels": np.zeros(4), "wheels_rpm": np.zeros(4),
            "magnetorquers": np.zeros(3), "rcs": np.zeros(8, dtype=bool)}


def saturate_commands(commands):
    """
    Clamps the ESP32 commands to physical limits:
      "ion_thruster" : float in [0, 1]   max thrust fraction (below min_hall_fraction -> off)
      "wheels"       : 4 floats in [-1,1] max torque fraction of each wheel
      "magnetorquers": 3 floats in [-1,1] max dipole fraction along body x, y, z
      "rcs"          : 8 booleans (on/off) in the order of corners (sx, sy, sz) = ++ +, ++ -, +- +, +- -, -+ +, -+ -, -- +, -- -
    """
    # Start with a safe baseline, then apply whatever the ESP32 sent.
    base = null_commands()
    base.update(commands)
    return {
        "ion_thruster": float(np.clip(base["ion_thruster"], 0.0, 1.0)),
        "wheels": np.clip(np.asarray(base["wheels"], dtype=float), -1.0, 1.0),
        "wheels_rpm": np.clip(np.asarray(base["wheels_rpm"], dtype=float),
                              -max_wheel_speed * 60.0 / (2.0 * np.pi), max_wheel_speed * 60.0 / (2.0 * np.pi)),
        "magnetorquers": np.clip(np.asarray(base["magnetorquers"], dtype=float), -1.0, 1.0),
        "rcs": np.asarray(base["rcs"], dtype=bool)
    }


def calculate_actuators(wheel_momentum, m_xenon, m_chemical, commands, R_b2i, B_eci):
    """
    Calculates forces, torques, and consumptions of the actuators in the body frame.
    Returns a dictionary with:
      body_force [N]            total thruster force (ion + RCS)
      body_torque [N m]         external torque on the satellite (thrusters + magnetorquers)
      wheel_torque [N m]        (4) torque applied to each wheel (after speed loop saturation)
      wheel_dh [N m]            (4) angular momentum derivative of the wheels
      dm_xenon, dm_chemical [kg/s]
    """
    # --- Hall thruster: variable intensity with minimum threshold ---
    u = commands["ion_thruster"]
    if u >= min_hall_fraction and m_xenon > 0.0:
        T = max_hall_thrust * u
    else:
        T = 0.0
    f_ion = T * HALL_DIRECTION
    mom_ion = np.cross(hall_position - center_of_mass_offset, f_ion)
    dm_xenon = -T / (g0 * isp_hall)

    # --- RCS: on/off, only if there is propellant ---
    # (If the tank is dry, these buttons are useless)
    active = commands["rcs"] if m_chemical > 0.0 else np.zeros(8, dtype=bool)
    f_rcs = UNIT_RCS_FORCE @ active.astype(float)
    mom_rcs = UNIT_RCS_TORQUE @ active.astype(float)
    dm_chemical = -np.sum(active) * rcs_thrust / (g0 * isp_rcs)

    # --- Magnetorquers: tau = m x B (in the body frame) ---
    dipole = commands["magnetorquers"] * max_mtq_dipole        # [A m^2]
    B_body = R_b2i.T @ B_eci                                   # [T]
    mom_mtq = np.cross(dipole, B_body)

    # --- Reaction wheels: torque saturates if the wheel is at its momentum limit ---
    # In SPEED mode, the wheel driver tracks the setpoint with a simple proportional loop,
    # limited by the max torque (the wheel can't accelerate faster than this).
    if WHEEL_MODE == "SPEED":
        omega_wheels = wheel_momentum / wheel_inertia                            # [rad/s]
        omega_cmd = commands["wheels_rpm"] * 2.0 * np.pi / 60.0                  # [rad/s]
        tau_w = np.clip(wheel_loop_gain * (omega_cmd - omega_wheels), -max_wheel_torque, max_wheel_torque)
    else:
        tau_w = commands["wheels"] * max_wheel_torque
        
    saturated = ((wheel_momentum >= max_wheel_momentum) & (tau_w > 0.0)) | ((wheel_momentum <= -max_wheel_momentum) & (tau_w < 0.0))
    tau_w = np.where(saturated, 0.0, tau_w)

    return {
        "body_force": f_ion + f_rcs,
        "body_torque": mom_ion + mom_rcs + mom_mtq,
        "wheel_torque": tau_w,
        "wheel_dh": tau_w,
        "dm_xenon": dm_xenon,
        "dm_chemical": dm_chemical
    }


# =====================================================================
# TOTAL ACCELERATION (All forces combined)
# =====================================================================
def calculate_total_acceleration(r_vec, v_vec, q, eph, m_xenon, m_chemical, actuator_body_force):
    """
    Calculates the total vector acceleration acting on the satellite,
    summing all physical contributions modeled so far, including actuator thrust.
    eph: dictionary returned by calculate_body_positions for the current instant.
    actuator_body_force: force [N] of Hall thruster + RCS in the body frame.
    """
    r_sun = eph["Sun"]
    r_moon = eph["Moon"]
    total_mass = dry_mass + m_xenon + m_chemical   # Instantaneous mass (propellant is consumed)

    # 1. Central Earth gravity (pure Keplerian)
    r_mag = np.linalg.norm(r_vec)
    a_total = - (mu / (r_mag ** 3)) * r_vec

    # Current attitude (needed for drag, SRP, and thrust)
    R_b2i = attitude_matrix(r_vec, v_vec, q)

    # 2. J2 Perturbation
    if PERTURBATIONS["J2"]:
        a_total = a_total + calculate_j2_acceleration(r_vec, mu, earth_equatorial_radius, J2)

    # 3. Third Body Perturbation - Sun
    if PERTURBATIONS["SUN"]:
        a_total = a_total + calculate_third_body_acceleration(r_vec, r_sun, mu_sun)

    # 4. Third Body Perturbation - Moon
    if PERTURBATIONS["MOON"]:
        a_total = a_total + calculate_third_body_acceleration(r_vec, r_moon, mu_moon)

    # 5. Atmospheric Drag
    if PERTURBATIONS["DRAG"]:
        a_total = a_total + calculate_drag_acceleration(
            r_vec, v_vec, total_mass,
            length, width, height,
            solar_panel_length, solar_panel_width,
            earth_equatorial_radius, earth_rotation_rate,
            R_b2i, r_sun
        )

    # 6. Solar Radiation Pressure (SRP)
    if PERTURBATIONS["SRP"]:
        a_total = a_total + calculate_srp_acceleration(
            r_vec, r_sun, total_mass,
            length, width, height,
            solar_panel_length, solar_panel_width,
            R_b2i
        )

    # 7. Actuator thrust (Hall thruster + RCS): N -> m/s^2 -> km/s^2, rotated to ECI
    a_total = a_total + (R_b2i @ actuator_body_force) / total_mass / 1000.0

    return a_total


# =====================================================================
# ATTITUDE DYNAMICS (used only in "FREE" mode)
# =====================================================================
def attitude_derivative(r_vec, q, w, wheel_momentum, m_total, external_body_torque, wheel_torque):
    """
    Quaternion and angular velocity derivatives (body frame) with reaction wheels.
    Total angular momentum: H = I w + A h_w. Equation:
        I dw/dt = tau_gg + tau_external - A * tau_wheels - w x (I w + A h_w)
    where tau_gg is the gravity gradient torque and A*tau_wheels the wheel reaction torque.
    """
    # Quaternion kinematics. Don't forget to normalize q elsewhere!
    q_dot = 0.5 * quat_multiply(q / np.linalg.norm(q), np.array([0.0, w[0], w[1], w[2]]))

    I = calculate_inertia_tensor(m_total)
    R_b2i = quat_to_matrix(q)
    r_mag = np.linalg.norm(r_vec)
    r_b = R_b2i.T @ (r_vec / r_mag)                                          # radial unit vector in the body frame
    torque_gg = 3.0 * mu / (r_mag ** 3) * np.cross(r_b, I @ r_b)             # [N m]

    H_wheels = WHEEL_AXES @ wheel_momentum                                   # wheel angular momentum [N m s]
    rhs = (torque_gg + external_body_torque
           - WHEEL_AXES @ wheel_torque
           - np.cross(w, I @ w + H_wheels))
    w_dot = np.linalg.solve(I, rhs)
    return q_dot, w_dot


# =====================================================================
# RUNGE-KUTTA 4 (RK4) NUMERICAL INTEGRATOR
# State y = [r (3), v (3), q (4), w (3), h_wheels (4), m_xenon, m_chemical]  -> 19 components
# =====================================================================
IDX_POS, IDX_VEL, IDX_QUAT, IDX_OMEGA, IDX_WHEELS, IDX_MXE, IDX_MCH = (slice(0, 3), slice(3, 6), slice(6, 10),
                                                                       slice(10, 13), slice(13, 17), 17, 18)


def state_derivative(y, eph, commands):
    r_vec, v_vec, q, w = y[IDX_POS], y[IDX_VEL], y[IDX_QUAT], y[IDX_OMEGA]
    wheel_momentum, m_xe, m_ch = y[IDX_WHEELS], y[IDX_MXE], y[IDX_MCH]

    R_b2i = attitude_matrix(r_vec, v_vec, q)
    B_eci = geomagnetic_field_eci(r_vec, eph["JD"])
    att = calculate_actuators(wheel_momentum, m_xe, m_ch, commands, R_b2i, B_eci)

    a = calculate_total_acceleration(r_vec, v_vec, q, eph, m_xe, m_ch, att["body_force"])

    if ATTITUDE_MODE == "FREE":
        m_tot = dry_mass + m_xe + m_ch
        q_dot, w_dot = attitude_derivative(r_vec, q, w, wheel_momentum, m_tot, att["body_torque"], att["wheel_torque"])
        h_dot = att["wheel_dh"]
    else:
        # Prescribed attitude mode: no dynamics to integrate here
        q_dot, w_dot, h_dot = np.zeros(4), np.zeros(3), np.zeros(4)

    return np.concatenate([v_vec, a, q_dot, w_dot, h_dot, [att["dm_xenon"], att["dm_chemical"]]])


def rk4_step(y, dt_seconds, current_dt, commands=None):
    """
    Performs a single RK4 integration step on the full state
    (position, velocity, attitude, wheels, masses) and returns the new state and time.
    Commands are kept constant during the step. Sun and Moon ephemerides are
    calculated 3 times per step (t, t+dt/2, t+dt).
    """
    commands = null_commands() if commands is None else saturate_commands(commands)
    dt_mid = current_dt + timedelta(seconds=dt_seconds / 2.0)
    dt_end = current_dt + timedelta(seconds=dt_seconds)

    # Calculate ephemerides at the three time points. It's a bit of an overhead
    # but ensures we don't miss fast changes in the Sun/Moon positions.
    eph0 = calculate_body_positions(current_dt, y[IDX_POS])
    eph1 = calculate_body_positions(dt_mid, y[IDX_POS])
    eph2 = calculate_body_positions(dt_end, y[IDX_POS])

    # Classic RK4 stages
    k1 = state_derivative(y, eph0, commands)
    k2 = state_derivative(y + 0.5 * dt_seconds * k1, eph1, commands)
    k3 = state_derivative(y + 0.5 * dt_seconds * k2, eph1, commands)
    k4 = state_derivative(y + dt_seconds * k3, eph2, commands)

    # Weighted average of the stages
    y_new = y + (dt_seconds / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    
    # Post-step corrections
    y_new[IDX_QUAT] /= np.linalg.norm(y_new[IDX_QUAT])                 # re-normalize the quaternion
    y_new[IDX_MXE] = max(y_new[IDX_MXE], 0.0)                          # tanks can't go negative
    y_new[IDX_MCH] = max(y_new[IDX_MCH], 0.0)
    
    return y_new, dt_end


def osculating_semi_major_axis(r_vec, v_vec):
    """Osculating semi-major axis [km] from specific energy."""
    energy = 0.5 * np.dot(v_vec, v_vec) - mu / np.linalg.norm(r_vec)
    return -mu / (2.0 * energy)


def eci_to_latlon(r_eci, jd):
    """Converts ECI [km] to (lon, lat) geographic [deg]. Spherical approximation (~10 km error)."""
    # GMST calculation for Earth rotation
    gmst = np.radians((280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360.0)
    cg, sg = np.cos(gmst), np.sin(gmst)
    
    # ECI -> ECEF rotation
    x =  cg * r_eci[0] + sg * r_eci[1]
    y = -sg * r_eci[0] + cg * r_eci[1]
    z =  r_eci[2]
    
    return float(np.degrees(np.arctan2(y, x))), float(np.degrees(np.arctan2(z, np.hypot(x, y))))


def nadir_pointing_error(y):
    """Angle [deg] between body +z axis (camera) and nadir direction."""
    R = attitude_matrix(y[IDX_POS], y[IDX_VEL], y[IDX_QUAT])
    nadir = -y[IDX_POS] / np.linalg.norm(y[IDX_POS])
    return np.degrees(np.arccos(np.clip(np.dot(R[:, 2], nadir), -1.0, 1.0)))


def initial_state():
    """Builds the full initial state from the configuration."""
    if USE_KEPLERIAN_ELEMENTS:
        r0, v0 = keplerian_to_eci(semi_major_axis, eccentricity, inclination, raan,
                                  argument_of_perigee, true_anomaly, mu)
    else:
        r0 = np.array([initial_pos_x, initial_pos_y, initial_pos_z], dtype=float)
        v0 = np.array([initial_vel_x, initial_vel_y, initial_vel_z], dtype=float)

    # Initial attitude: align with LVLH (nadir pointing) plus any configured offset
    q0 = matrix_to_quat(lvlh_matrix(r0, v0) @ ATTITUDE_OFFSET_MATRIX)
    
    w0 = np.array([initial_omega_x, initial_omega_y, initial_omega_z], dtype=float)
    if INITIAL_OMEGA_NADIR:
        # If we want to start nadir-pointing, we must add the orbital angular velocity
        n_orb = np.linalg.norm(np.cross(r0, v0)) / np.linalg.norm(r0) ** 2        # orbital angular velocity [rad/s]
        w0 = w0 + ATTITUDE_OFFSET_MATRIX.T @ np.array([0.0, -n_orb, 0.0])
        
    return np.concatenate([r0, v0, q0, w0, np.asarray(initial_wheel_momentum, dtype=float),
                           [xenon_mass, chemical_propellant_mass]])


# =====================================================================
# HWIL COMMUNICATION WITH THE ESP32 (line-based JSON) AND SIMULATION LOOP
# The PC simulates ONLY the physics and sensors. The decisions (when to fire an actuator)
# are made by the ESP32.
# =====================================================================

# --- Parameters (CUSTOMIZABLE) ---
SERIAL_PORT = "/dev/cu.usbserial-110"   # Change this to "COM5" on Windows, or "/dev/ttyUSB0" on Linux
BAUD_RATE = 115200                      # Must match Serial.begin() in the firmware
dt_sim = 1.0                            # [s] Simulated step/tick (sensor periods must be multiples of this)
k_factor = 25.0                         # [-] Simulated seconds per real second
max_real_time = 365 * 24 * 60 * 60      # [s] Max real duration (effectively no limit for now)
ESP32_TIMEOUT_S = None                  # [sim s] No messages from ESP32 for this long -> safe mode (None = disabled)
max_simulated_time = None               # [s] Max simulated duration (None = no limit)
logging_interval_s = 60.0               # [sim s] Recording interval for the history

# --- ASCII dashboard ---
LIVE_DASHBOARD = True                   # True: draw the dashboard during the simulation
DASHBOARD_UPDATE_INTERVAL_S = 1.0       # [real s] Dashboard refresh rate
                                        # (independent of k; 1 line per real second)
                                        # Note: with LIVE_DASHBOARD=True, the textual ESP32 logs
                                        # get overwritten by the dashboard on every refresh.

# --- Live map (in a separate process) ---
MAP_ACTIVE = True                       # True: native window, map updated
MAP_UPDATE_INTERVAL_S = 2.0             # [real s] Map update rate
MAP_MAX_POINTS = 500                    # Ground track points kept (~9 orbits)


# Sensor frequencies (simulated time).
# For long simulations, set dt_sim = 10 and keep only "gps" and "tanks".
SENSOR_FREQUENCIES_HZ = {
    "gps": 0.1,              # -> "pos", "vel"
    "gyroscope": 1.0,        # -> "gyro"
    "star_tracker": 0.5,     # -> "q"
    "magnetometer": 1.0,     # -> "mag"
    "sun_sensor": 0.5,       # -> "sun", "eclipse"
    "wheels": 1.0,           # -> "rw"
    "tanks": 0.1,            # -> "tank"
}


# --- Sensors: true state -> JSON fields (NO noise: you can add it later) ---
def _fmt(v, sig_figs):
    """Rounds to 'sig_figs' significant digits (reduces bytes sent)."""
    return [float(f"{x:.{sig_figs}g}") for x in np.atleast_1d(v)]


def true_attitude(y):
    """(Body->ECI R matrix, quaternion [w,x,y,z], body angular velocity) true values, valid in both modes."""
    r, v = y[IDX_POS], y[IDX_VEL]
    R = attitude_matrix(r, v, y[IDX_QUAT])
    if ATTITUDE_MODE == "FREE":
        return R, y[IDX_QUAT] / np.linalg.norm(y[IDX_QUAT]), y[IDX_OMEGA]
    
    # Prescribed LVLH attitude: angular velocity = orbital rate
    n = np.linalg.norm(np.cross(r, v)) / np.dot(r, r)
    return R, matrix_to_quat(R), ATTITUDE_OFFSET_MATRIX.T @ np.array([0.0, -n, 0.0])


def sample_sensor(name, y, eph):
    R, q, w = true_attitude(y)
    
    if name == "gyroscope":
        return {"gyro": _fmt(w, 6)}                                                   # [rad/s]
        
    if name == "magnetometer":
        return {"mag": _fmt(R.T @ geomagnetic_field_eci(y[IDX_POS], eph["JD"]), 5)}   # [T], body frame
        
    if name == "star_tracker":
        return {"q": _fmt(q, 8)}                                                      # [w,x,y,z] body -> ECI
        
    if name == "sun_sensor":
        s = eph["Sun"] - y[IDX_POS]
        return {"sun": _fmt(R.T @ (s / np.linalg.norm(s)), 6),                        # body frame unit vector
                "eclipse": int(illumination_factor(y[IDX_POS], eph["Sun"]) < 0.5)}    # 0 = sun, 1 = eclipse
                
    if name == "gps":
        return {"pos": _fmt(y[IDX_POS], 9), "vel": _fmt(y[IDX_VEL], 9)}               # ECI [km], [km/s]
        
    if name == "wheels":
        return {"rw": [round(float(x), 2) for x in y[IDX_WHEELS] / wheel_inertia * 60.0 / (2.0 * np.pi)]}   # [RPM]
        
    if name == "tanks":
        return {"tank": [round(float(y[IDX_MXE]), 4), round(float(y[IDX_MCH]), 4)]}   # [kg]
        
    raise ValueError(f"Unknown sensor: {name}")


def generate_packet(y, current_dt, t_sim, tick, tick_periods):
    """A single JSON with "t" and ONLY the sensors sampling in this tick (None if none)."""
    due_sensors = [name for name, period in tick_periods.items() if tick % period == 0]
    if not due_sensors:
        return None
        
    eph = calculate_body_positions(current_dt, y[IDX_POS])
    packet = {"t": round(t_sim, 3)}
    for name in due_sensors:
        packet.update(sample_sensor(name, y, eph))
        
    return packet


# --- ESP32 Commands -> Actuators (PERSISTENT state: stays as is until a change arrives) ---
def initial_commands(y):
    """Initial condition: everything off; wheels maintain their starting speed."""
    c = null_commands()
    c["wheels_rpm"] = y[IDX_WHEELS] / wheel_inertia * 60.0 / (2.0 * np.pi)
    return c


def apply_esp32_message(msg, commands):
    """Updates commands with ONLY the present and VALID fields. Returns the number of discarded fields."""
    discarded = 0

    def parse_vector(val, n):
        a = np.asarray(val, dtype=float)
        if a.shape != (n,) or not np.all(np.isfinite(a)):
            raise ValueError("invalid shape or values")
        return a

    if "hall" in msg:
        try:
            h = float(msg["hall"])
            if not np.isfinite(h):
                raise ValueError("hall is not finite")
            commands["ion_thruster"] = h / 100.0                                 # [%] -> [0,1]
        except (TypeError, ValueError):
            discarded += 1
            
    if "rw" in msg:
        try:
            commands["wheels_rpm"] = parse_vector(msg["rw"], 4)                  # [RPM] with sign
        except (TypeError, ValueError):
            discarded += 1
            
    if "mtq" in msg:
        try:
            commands["magnetorquers"] = parse_vector(msg["mtq"], 3) / max_mtq_dipole   # [A m^2] -> [-1,1]
        except (TypeError, ValueError):
            discarded += 1
            
    if "rcs" in msg:
        try:
            commands["rcs"] = parse_vector(msg["rcs"], 8) != 0.0                 # list of 8 0/1 values
        except (TypeError, ValueError):
            discarded += 1
            
    return discarded


# --- Serial connection (non-blocking read) ---
class SerialConnection:
    MAX_BUFFER_SIZE = 4096               # [bytes] above this, without "\n", the buffer is cleared

    def __init__(self, port=SERIAL_PORT, baud=BAUD_RATE):
        import serial
        self._write_timeout_exc = serial.SerialTimeoutException
        self.ser = serial.Serial()
        self.ser.port, self.ser.baudrate, self.ser.timeout, self.ser.write_timeout = port, baud, 0, 1
        
        # Toggling DTR/RTS avoids resetting the ESP32 when opening the port.
        # (On some boards, it resets anyway, which is fine but delays the boot).
        self.ser.dtr = False
        self.ser.rts = False
        self.ser.open()
        
        time.sleep(2.0)                  # Wait for the bootloader to finish
        self.ser.reset_input_buffer()
        
        self._buf = b""
        self.bytes_sent = 0
        self.bytes_received = 0
        self.skipped_lines = 0           # Start with "{" but aren't valid JSON objects
        self.ignored_lines = 0           # Text lines (boot logs, debug prints)
        self.write_timeouts = 0

    def send(self, dictionary):
        """Sends a JSON line. True if written, False if the write timed out."""
        line = (json.dumps(dictionary, separators=(",", ":")) + "\n").encode()
        try:
            self.ser.write(line)
        except self._write_timeout_exc:
            self.write_timeouts += 1
            self.ser.reset_output_buffer()
            return False
        
        self.bytes_sent += len(line)
        return True

    def read_messages(self):
        """List of (arrival_time, message) received so far. Never blocks."""
        n = self.ser.in_waiting
        if n:
            data = self.ser.read(n)
            self.bytes_received += len(data)
            self._buf += data
            
        arrival_time = time.perf_counter()
        
        # Buffer overflow protection: if we get a huge chunk of garbage, drop it.
        if len(self._buf) > self.MAX_BUFFER_SIZE and b"\n" not in self._buf:
            self._buf = b""
            self.skipped_lines += 1
            
        messages = []
        while b"\n" in self._buf:
            line, self._buf = self._buf.split(b"\n", 1)
            line = line.strip()
            if not line:
                continue
                
            # Ignore anything that doesn't look like a JSON object
            if not line.startswith(b"{"):
                self.ignored_lines += 1
                continue
                
            try:
                msg = json.loads(line)
            except ValueError:
                self.skipped_lines += 1
                continue
                
            if isinstance(msg, dict):
                messages.append((arrival_time, msg))
            else:
                self.skipped_lines += 1
                
        return messages

    def close(self):
        self.ser.close()


def print_dashboard(y, commands, t_sim, tick, current_dt,
                    n_msg_esp, latencies_ms, safe_mode,
                    late_cycles, real_elapsed_time,
                    panel_integrity=1.0, total_impacts=0,
                    power_proxy_W=SOLAR_PANEL_POWER_BOL_W, expected_lambda_tick=0.0,
                    conjunctions=None, n_cam_executed=0, n_cam_successful=0, n_collisions=0, ping_status="waiting", n_pings_received=0):
    if conjunctions is None:
        conjunctions = []
    
    """Full-screen ASCII dashboard: called from the run_simulation loop."""
    r_vec, v_vec = y[IDX_POS], y[IDX_VEL]

    # --- Osculating orbit ---
    r_mag = float(np.linalg.norm(r_vec))
    v2    = float(np.dot(v_vec, v_vec))
    a_o   = -mu / (2.0 * (0.5 * v2 - mu / r_mag))
    
    # Moving average over 1 orbit (perfectly filters out J2)
    if not hasattr(print_dashboard, "_a_buf"):
        print_dashboard._a_buf = []
    print_dashboard._a_buf.append((t_sim, a_o))
    sim_period = 2.0 * np.pi * np.sqrt(a_o**3 / mu)
    while print_dashboard._a_buf and print_dashboard._a_buf[0][0] < t_sim - sim_period:
        print_dashboard._a_buf.pop(0)
    a_mean = float(np.mean([v for _, v in print_dashboard._a_buf]))
    
    # Only show the mean if the buffer covers ~1 period (otherwise it's biased by J2)
    buf_duration = t_sim - print_dashboard._a_buf[0][0] if print_dashboard._a_buf else 0.0
    if buf_duration < 0.99 * sim_period:
        a_mean = float('nan')   # will be printed as "nan" -> you know it's a transient
        
    e_vec = ((v2 - mu / r_mag) * r_vec - np.dot(r_vec, v_vec) * v_vec) / mu
    ecc   = float(np.linalg.norm(e_vec))
    perigee = a_o * (1.0 - ecc)
    apogee  = a_o * (1.0 + ecc)
    period = 2.0 * np.pi * np.sqrt(a_o ** 3 / mu) / 60.0      # min
    altitude = r_mag - earth_equatorial_radius                # km (spherical, display only)

    # --- Attitude ---
    pointing_err = nadir_pointing_error(y) if ATTITUDE_MODE == "FREE" else 0.0
    w = y[IDX_OMEGA]

    # --- Illumination (requires ephemerides) ---
    eph = calculate_body_positions(current_dt, r_vec)
    illum = illumination_factor(r_vec, eph["Sun"])
    if   illum >= 0.99:  sun_state = "FULL SUN "
    elif illum >= 0.50:  sun_state = "PENUMBRA "
    else:                sun_state = "ECLIPSE  "

    # --- Tanks ---
    m_xe, m_ch = float(y[IDX_MXE]), float(y[IDX_MCH])
    pct_xe = 100.0 * m_xe / xenon_mass if xenon_mass > 0 else 0.0
    pct_ch = 100.0 * m_ch / chemical_propellant_mass if chemical_propellant_mass > 0 else 0.0

    # --- Wheels (actual and commanded, RPM) ---
    rpm_actual = y[IDX_WHEELS] / wheel_inertia * 60.0 / (2.0 * np.pi)
    rpm_cmd    = commands["wheels_rpm"]

    # --- Actuators ---
    # --- HALL: percentage + effective thrust + mass flow ---
    hall_pct = commands["ion_thruster"] * 100.0
    if commands["ion_thruster"] >= min_hall_fraction and y[IDX_MXE] > 0.0:
        T_hall = max_hall_thrust * commands["ion_thruster"]
        dm_xe  = T_hall / (g0 * isp_hall)
    else:
        T_hall, dm_xe = 0.0, 0.0
    hall_str = (f"HALL {hall_pct:5.1f}%  (thrust {T_hall*1e3:5.2f} mN, "
                f"dm {dm_xe:7.2e} kg/s)") if hall_pct > 0.5 else "HALL OFF"

    # --- MTQ: 3 components + magnitude ---
    mtq_x, mtq_y, mtq_z = (commands["magnetorquers"] * max_mtq_dipole)
    mtq_norm = float(np.sqrt(mtq_x**2 + mtq_y**2 + mtq_z**2))
    mtq_str = (f"MTQ [A·m²]: m = [{mtq_x:+6.2f} {mtq_y:+6.2f} {mtq_z:+6.2f}]  "
               f"|m| = {mtq_norm:5.2f}")

    # --- RCS: bitmask + which ones are active ---
    rcs_bits = commands["rcs"]                          # array of 8 bools
    rcs_visual = "".join("#" if b else "." for b in rcs_bits)
    rcs_mask = sum(int(b) << i for i, b in enumerate(rcs_bits))
    n_rcs = int(np.sum(rcs_bits))
    rcs_str = (f"RCS: {rcs_visual}  (0x{rcs_mask:02X}, {n_rcs}/8)  "
               f"order: +++ ++- +-+ +-- -++ -+- --+ ---")

    # --- Wheels: saturation h/h_max (%) ---
    wheel_sat = 100.0 * y[IDX_WHEELS] / max_wheel_momentum    # [-100%, +100%]

    # --- ESP32 Link ---
    if len(latencies_ms) > 0:
        arr = np.asarray(latencies_ms)
        lat_str = f"med {np.median(arr):5.1f} ms  p95 {np.percentile(arr, 95):5.1f} ms"
    else:
        lat_str = "no measurements (yet)"

    # --- Percentage bars ---
    def draw_bar(pct, width=24):
        n = int(round(max(0.0, min(100.0, pct)) / 100.0 * width))
        return "#" * n + "." * (width - n)

    # --- Clear screen in Jupyter (no-op from terminal) ---
    clear_output(wait=True)

    # --- Print (I know, this formatting is a pain, but it looks good on a terminal) ---
    print("=" * 74)
    print("        SATELLITE SIMULATOR — LIVE DASHBOARD (PC + ESP32)        ")
    print("=" * 74)
    print(f" SIM TIME   {t_sim/3600:8.3f} h ({t_sim:7.0f} s)   |  TICK {tick:7d}")
    print(f" T-ZERO     {current_dt.strftime('%Y-%m-%d %H:%M:%S')}")
    if real_elapsed_time > 0:
        print(f" REAL       {real_elapsed_time:6.1f} s elapsed   |  "
              f"effective k {t_sim/real_elapsed_time:5.1f}")
    else:
        print(f" REAL       0.0 s elapsed")
    print("-" * 74)
    print(" [ORBIT]")
    print(f"   a (osc) = {a_o:10.3f} km    e (osc) = {ecc:8.5f}")
    if np.isnan(a_mean):
        print(f"   a (mean) =   in transient (first orbit)")
    else:
        print(f"   a (mean) = {a_mean:8.3f} km   Δ vs nominal = {a_mean - 6878.14:+7.3f} km")
    print(f"   perigee = {perigee:10.3f} km    apogee = {apogee:10.3f} km")
    print(f"   period  = {period:10.2f} min   alt~ = {altitude:10.1f} km")
    print("-" * 74)
    print(" [ATTITUDE]")
    print(f"   Nadir pointing error : {pointing_err:8.4f} deg")
    print(f"   Illumination         : {sun_state}  (nu = {illum:5.3f})")
    print(f"   body omega [rad/s]   : [{w[0]:+9.5f} {w[1]:+9.5f} {w[2]:+9.5f}]")
    print("-" * 74)
    print(" [CONSUMABLES]")
    print(f"   Xenon   : {pct_xe:5.1f}%  [{draw_bar(pct_xe)}]  {m_xe:7.4f} kg")
    print(f"   Chem    : {pct_ch:5.1f}%  [{draw_bar(pct_ch)}]  {m_ch:7.4f} kg")
    print("-" * 74)
    print("-" * 74)
    print(" [DEBRIS / PANELS]")
    pct_int = 100.0 * panel_integrity
    print(f"   Panel integrity : {pct_int:6.2f}%  [{draw_bar(pct_int)}]")
    print(f"   Electrical power: {power_proxy_W:7.1f} W   "
          f"(BOL: {SOLAR_PANEL_POWER_BOL_W:.0f} W, loss: "
          f"{SOLAR_PANEL_POWER_BOL_W - power_proxy_W:.2f} W)")
    print(f"   Cumulative hits : {total_impacts:6d}   "
          f"(λ_physical = {expected_lambda_tick:.3e} impacts/tick)")
    print("-" * 74)
    print(" [LARGE DEBRIS]")
    print(f"   Active conjunctions : {len(conjunctions)}")
    print(f"   CAM executed/success: {n_cam_executed}/{n_cam_successful}")
    print(f"   Collisions          : {n_collisions}")
    for c in conjunctions:
        if c["state"] == "in_progress":
            print(f"     -> id={c['id']} TCA in {(c['tca']-t_sim)/60:.0f} min, r_target {c['r_target_km']:.2f} km")
    print(" [ACTUATORS]")
    print(f"   {hall_str}")
    print(f"   {mtq_str}")
    print(f"   {rcs_str}")
    print(f"   Wheels RPM cmd  : [{rpm_cmd[0]:+8.1f} {rpm_cmd[1]:+8.1f} "
          f"{rpm_cmd[2]:+8.1f} {rpm_cmd[3]:+8.1f}]")
    print(f"   Wheels RPM act  : [{rpm_actual[0]:+8.1f} {rpm_actual[1]:+8.1f} "
          f"{rpm_actual[2]:+8.1f} {rpm_actual[3]:+8.1f}]")
    print(f"   Wheels sat [%]  : [{wheel_sat[0]:+8.1f} {wheel_sat[1]:+8.1f} "
          f"{wheel_sat[2]:+8.1f} {wheel_sat[3]:+8.1f}]")
    print("-" * 74)
    print(" [ESP32 LINK]")
    print(f"   Messages received : {n_msg_esp}")
    print(f"   Ping/pong         : {ping_status}   ({n_pings_received} total replies)")
    print(f"   RTT latency       : {lat_str}")
    print(f"   Late cycles       : {late_cycles}")
    if safe_mode:
        print("   *** SAFE MODE ACTIVE (ESP32 silent) ***")
    print("=" * 74)
    print(" (Interrupt the kernel in Jupyter to stop the simulation)")


_stop_requested = False

def _sigint_handler(signum, frame):
    """Ctrl+C: requests the loop to stop, so results aren't lost."""
    global _stop_requested
    _stop_requested = True


# --- Main Loop ---
def run_simulation(connection, k=None, max_real_time_arg=-1, max_sim_time_arg=-1, verbose=True):
    k = k_factor if k is None else k
    max_real_time_arg = max_real_time if max_real_time_arg == -1 else max_real_time_arg
    max_sim_time_arg = max_simulated_time if max_sim_time_arg == -1 else max_sim_time_arg

    # Sanity check: if the wheel loop is too fast for our simulation step, it will explode.
    if ATTITUDE_MODE == "FREE" and WHEEL_MODE == "SPEED" and dt_sim > wheel_loop_tau / 2.0:
        raise ValueError(f"dt_sim ({dt_sim} s) too large for the wheel loop (tau = {wheel_loop_tau} s): "
                         f"use dt_sim <= {wheel_loop_tau / 2.0} s or set ATTITUDE_MODE = 'LVLH'")

    y = initial_state()
    current_dt = datetime(year, month, day, hour, minute, second)
    commands = initial_commands(y)                                        # INITIAL CONDITION: all off

    tick_periods = {}
    for name, f in SENSOR_FREQUENCIES_HZ.items():
        n = (1.0 / f) / dt_sim
        if abs(n - round(n)) > 1e-9 or round(n) < 1:
            raise ValueError(f"Sensor '{name}' period ({1.0 / f} s) is not a multiple of dt_sim ({dt_sim} s)")
        tick_periods[name] = int(round(n))

    connection.send({"cmd": "reset"})                                    # Reset the on-board computer state
    time.sleep(0.2)
    connection.read_messages()

    cycle_duration = dt_sim / k                                          # Target real time per tick
    t_start = time.perf_counter()
    deadline = t_start                                                   # ABSOLUTE deadline for the current cycle
    t_sim, tick = 0.0, 0
    late_cycles, max_delay = 0, 0.0
    
    n_msg_esp, cmd_delays, latencies_ms, pointing_errors, events = 0, [], [], [], []
    send_times, msg_queue = {}, []                                       # sim time -> real time of sending
    
    # --- Periodic ESP32 link check via ping/pong (REAL time) ---
    PING_INTERVAL_S = 60.0          # How many REAL seconds between pings
    PING_TIMEOUT_S = 5.0            # Beyond this time without pong -> marked as TIMEOUT
    last_ping_real_time = -PING_INTERVAL_S   # force a ping on the first tick
    pending_ping = None             # real time of the ping awaiting pong (None = no open ping)
    ping_status = "waiting"         # text shown in the dashboard
    n_pings_received = 0            # cumulative pong replies
    
    last_msg_time, safe_mode, discarded_fields = 0.0, False, 0
    history = {"t": [], "pointing_error": [], "osc_perigee": [], "ion_thruster": [], "xenon": [],
               "wheel_rpm": [], "panel_integrity": [], "total_impacts": [],
               "power_proxy_W": [], "expected_lambda": [],
               "a_osc": [], "e_osc": [], "i_osc_deg": [], "LTAN_h": []}
    last_log_time = -1e18
    
    # --- Debris / Panels state (outside state vector y: it's a jump process) ---
    panel_integrity = INITIAL_PANEL_INTEGRITY
    total_impacts = 0
    impacts_tick = 0
    expected_lambda_tick = 0.0
    debris_rng = np.random.default_rng(DEBRIS_SEED)
    panel_area_m2 = 2.0 * solar_panel_length * solar_panel_width         # 1.8 m²
    power_proxy_W = SOLAR_PANEL_POWER_BOL_W * panel_integrity
    
    last_dashboard_time = -1e18             # [real s] last dashboard refresh
    
    # Reset the moving average buffer for 'a' (function attribute, persists across runs)
    if hasattr(print_dashboard, "_a_buf"):
        print_dashboard._a_buf = []
        
    # --- Live map in SEPARATE PROCESS (never blocks the main loop) ---
    _map_queue, _map_stop, _map_proc = None, None, None
    last_map_time = -1e18
    if MAP_ACTIVE and _HAS_MAP_PROC and _HAS_CARTOPY:
        try:
            _orbit_period_s = 2.0 * np.pi * np.sqrt(semi_major_axis ** 3 / mu)
            _map_queue = mp.Queue(maxsize=20)
            _map_stop = mp.Event()
            _map_proc = mp.Process(
                target=_map_process_fn,
                args=(_map_queue, MAP_MAX_POINTS, _map_stop, _orbit_period_s),
                daemon=True,
            )
            _map_proc.start()
            print(f"Live map started in separate process (PID {_map_proc.pid}).")
            print("Close the map window to close it; the simulation keeps running.")
        except Exception as e:
            print(f"Live map not available: {e}")
            _map_proc = None
    elif MAP_ACTIVE:
        print("⚠️  Map not active: cartopy or map_process.py is missing.")
        
    # --- Large debris: volume, effective area, conjunction queue ---
    satellite_volume_m3 = (width * (length + 2.0 * solar_panel_length) * height)
    effective_area_m2 = satellite_volume_m3 ** (2.0 / 3.0)
    large_debris_rng = np.random.default_rng(DEBRIS_SEED + 1 if DEBRIS_SEED is not None else None)
    conjunctions = []   # list of dicts: id, tca, start_radius, target_radius, state
    next_conjunction_id = 0
    n_cam_executed = 0
    n_cam_successful = 0
    n_collisions = 0


    while True:
        # Check if the user hit Ctrl+C or if we reached the time limits
        if _stop_requested:
            reason = "interrupted by user"; break
        if max_real_time_arg is not None and time.perf_counter() - t_start >= max_real_time_arg:
            reason = "max real time reached"; break
        if max_sim_time_arg is not None and t_sim >= max_sim_time_arg:
            reason = "max simulated time reached"; break

        # 1) Read incoming commands from ESP32 (NON-blocking; many might have arrived during the last wait)
        msg_queue.extend(connection.read_messages())
        for arrival_time, msg in msg_queue:
            n_msg_esp += 1
            last_msg_time = t_sim
            
            if safe_mode:
                safe_mode = False
                events.append((t_sim, "ESP32 is back online: exiting safe mode"))
                
            discarded_fields += apply_esp32_message(msg, commands)
            
            # If the message has a timestamp, calculate the latency
            if "t" in msg:
                try:
                    echo_t = round(float(msg["t"]), 3)
                except (TypeError, ValueError):
                    echo_t = None
                if echo_t is not None:
                    cmd_delays.append(t_sim - echo_t)  # delay in SIM seconds
                    send_time = send_times.pop(echo_t, None)
                    if send_time is not None:
                        latencies_ms.append((arrival_time - send_time) * 1e3)  # round-trip in REAL ms
                        
            # Handle CAM (Collision Avoidance Maneuver) responses
            if "cam" in msg:
                cam = msg["cam"]
                cid = cam.get("id")
                cam_state = cam.get("status")
                for c in conjunctions:
                    if c["id"] == cid:
                        if cam_state == "complete":
                            c["state"] = "avoided"
                            n_cam_executed += 1
                            n_cam_successful += 1
                            events.append((t_sim, f"CAM id={cid} completed (altitude reached)"))
                        elif cam_state == "aborted":
                            events.append((t_sim, f"CAM id={cid} aborted"))
                        break
                        
            # Handle ping responses
            if msg.get("log") == "pong" and pending_ping is not None:
                rtt_ms = (arrival_time - pending_ping) * 1e3
                ping_status = f"OK (RTT {rtt_ms:.1f} ms)"
                n_pings_received += 1
                pending_ping = None
                
            # Log event
            if "log" in msg:
                events.append((t_sim, msg["log"]))
                if verbose:
                    print(f"  [ESP32 @ t={t_sim / 3600:.3f} h] {msg['log']}")
        msg_queue.clear()

        # 1c) Ping timeout check
        if pending_ping is not None and (time.perf_counter() - pending_ping) > PING_TIMEOUT_S:
            ping_status = f"TIMEOUT (>{PING_TIMEOUT_S:.0f} s without pong)"
            pending_ping = None

        # 1b) Watchdog: if the ESP32 goes silent for too long, we assume something broke.
        if ESP32_TIMEOUT_S is not None and not safe_mode and t_sim - last_msg_time > ESP32_TIMEOUT_S:
            safe_mode = True
            commands["ion_thruster"] = 0.0
            commands["rcs"] = np.zeros(8, dtype=bool)
            commands["magnetorquers"] = np.zeros(3)
            commands["wheels_rpm"] = y[IDX_WHEELS] / wheel_inertia * 60.0 / (2.0 * np.pi)   # wheels hold their speed
            events.append((t_sim, "WATCHDOG: no messages from ESP32, actuators in safe mode"))
            if verbose:
                print(f"  [WATCHDOG @ t={t_sim / 3600:.3f} h] ESP32 silent: actuators in safe mode")

        # 1d) Periodic ping (REAL time, every PING_INTERVAL_S seconds)
        # Just to check if the ESP32 is alive and the serial link works,
        # regardless of the actuators (which usually stay silent if nothing changes).
        if time.perf_counter() - last_ping_real_time >= PING_INTERVAL_S:
            last_ping_real_time = time.perf_counter()
            if pending_ping is not None:
                ping_status = "TIMEOUT (previous ping unanswered)"
            connection.send({"cmd": "ping"})
            pending_ping = last_ping_real_time

        # 2) Sensors -> ESP32 (only those sampling in this tick, each at its own frequency)
        sensor_packet = generate_packet(y, current_dt, t_sim, tick, tick_periods)
        if sensor_packet:
            if connection.send(sensor_packet):
                send_times[sensor_packet["t"]] = time.perf_counter()
                if len(send_times) > 1000:  # prevent unbounded growth if the ESP32 is dead
                    send_times.pop(next(iter(send_times)))
                    
        # 2a-bis) Live map: send position via queue (never blocks)
        if _map_proc is not None and _map_proc.is_alive() and \
           (time.perf_counter() - last_map_time) >= MAP_UPDATE_INTERVAL_S:
            last_map_time = time.perf_counter()
            try:
                _eph_map = calculate_body_positions(current_dt, y[IDX_POS])
                jd_map = _eph_map["JD"]
                lon_s, lat_s = eci_to_latlon(y[IDX_POS], jd_map)
                
                # Subsolar point = Sun direction converted to geographic lat/lon
                lon_sun, lat_sun = eci_to_latlon(_eph_map["Sun"], jd_map)
                try:
                    # (t_sim, lon_sat, lat_sat, sub_lat, sub_lon)
                    _map_queue.put_nowait((t_sim, lon_s, lat_s, lat_sun, lon_sun))
                except Exception:
                    pass    # queue full (child is slow): skip this frame
            except Exception:
                pass

        # 2b) Live ASCII dashboard (refresh every DASHBOARD_UPDATE_INTERVAL_S REAL seconds)
        if LIVE_DASHBOARD and (time.perf_counter() - last_dashboard_time) >= DASHBOARD_UPDATE_INTERVAL_S:
            last_dashboard_time = time.perf_counter()
            print_dashboard(y, commands, t_sim, tick, current_dt,
                            n_msg_esp, latencies_ms, safe_mode,
                            late_cycles, time.perf_counter() - t_start,
                            panel_integrity=panel_integrity,
                            total_impacts=total_impacts,
                            power_proxy_W=power_proxy_W,
                            expected_lambda_tick=expected_lambda_tick,
                            conjunctions=conjunctions,
                            n_cam_executed=n_cam_executed,
                            n_cam_successful=n_cam_successful,
                            n_collisions=n_collisions,
                            ping_status=ping_status,
                            n_pings_received=n_pings_received)

        # 2b-bis) Large debris: warning generation and maneuver handling
        current_altitude_km = geodetic_altitude(y[IDX_POS])  # Calculated here, used below
        if LARGE_DEBRIS_ACTIVE:
            rho_large = get_large_debris_density(current_altitude_km)
            lam_large = (rho_large * (effective_area_m2 * 1e-6) *
                         LARGE_DEBRIS_REL_V_KM_S * dt_sim *
                         LARGE_DEBRIS_ACCELERATION)
            n_events = int(large_debris_rng.poisson(lam_large))
            
            for _ in range(n_events):
                tca = t_sim + CAM_WARNING_TIME_S
                start_radius_km = float(np.linalg.norm(y[IDX_POS]))
                c = {"id": next_conjunction_id,
                     "tca": tca,
                     "start_radius_km": start_radius_km,
                     "r_target_km": start_radius_km + CAM_ALTITUDE_DELTA_KM,
                     "state": "in_progress"}   # in_progress -> avoided / collision
                conjunctions.append(c)
                
                # Send the warning to the ESP32 (JSON "conj" message)
                conj_packet = {"t": round(t_sim, 3),
                               "conj": {"id": c["id"], "tca": round(tca, 1),
                                        "r_target_km": round(c["r_target_km"], 3)}}
                if connection.send(conj_packet):
                    send_times[conj_packet["t"]] = time.perf_counter()
                next_conjunction_id += 1

            # Pre-TCA check: mark as "avoided" as soon as |r| >= target
            for c in conjunctions[:]:
                if c["state"] == "in_progress":
                    current_radius = float(np.linalg.norm(y[IDX_POS]))
                    if current_radius >= c["r_target_km"]:
                        c["state"] = "avoided"
                        n_cam_executed += 1
                        n_cam_successful += 1
                        events.append((t_sim, f"CAM id={c['id']} AVOIDED (pre-TCA): "
                                             f"|r|={current_radius:.3f} >= {c['r_target_km']:.3f} km"))
                        conjunctions.remove(c)

            # Check if any conjunction reached the TCA
            for c in conjunctions[:]:
                if c["state"] == "in_progress" and t_sim >= c["tca"]:
                    current_radius = float(np.linalg.norm(y[IDX_POS]))
                    if current_radius >= c["r_target_km"]:
                        c["state"] = "avoided"
                        n_cam_executed += 1
                        n_cam_successful += 1
                        events.append((t_sim, f"CAM id={c['id']} AVOIDED: |r|={current_radius:.3f} >= {c['r_target_km']:.3f} km"))
                    elif current_radius >= c["r_target_km"] - 0.5:
                        # CAM failed but satellite is safe (within 500m of target)
                        c["state"] = "avoided"
                        n_cam_executed += 1
                        # not counted as successful, but it's not a collision
                        events.append((t_sim, f"CAM id={c['id']} FAILED (margin): |r|={current_radius:.3f} vs target {c['r_target_km']:.3f} km"))
                    else:
                        # Physical collision: way too far from the target
                        c["state"] = "collision"
                        n_collisions += 1
                        panel_integrity = 0.0
                        power_proxy_W = 0.0
                        events.append((t_sim, f"COLLISION id={c['id']}: |r|={current_radius:.3f} << {c['r_target_km']:.3f} km"))
                    conjunctions.remove(c)
                elif c["state"] == "avoided":
                    conjunctions.remove(c)

        # 2c) Small debris: update panel integrity (Poisson process) and power proxy
        panel_integrity, total_impacts, impacts_tick, expected_lambda_tick = \
            update_panel_integrity(
                panel_integrity, total_impacts,
                current_altitude_km, panel_area_m2,
                dt_sim, debris_rng)
        power_proxy_W = SOLAR_PANEL_POWER_BOL_W * panel_integrity

        # 3) Recording (true state, for analysis only)
        if ATTITUDE_MODE == "FREE":
            pointing_errors.append(nadir_pointing_error(y))
            
        if t_sim - last_log_time >= logging_interval_s:
            last_log_time = t_sim
            rn = np.linalg.norm(y[IDX_POS]); v2 = np.dot(y[IDX_VEL], y[IDX_VEL])
            a_o = -mu / (2.0 * (0.5 * v2 - mu / rn))
            e_vec = ((v2 - mu / rn) * y[IDX_POS] - np.dot(y[IDX_POS], y[IDX_VEL]) * y[IDX_VEL]) / mu
            
            history["t"].append(t_sim)
            history["pointing_error"].append(pointing_errors[-1] if pointing_errors else 0.0)
            history["osc_perigee"].append(a_o * (1.0 - np.linalg.norm(e_vec)))
            history["ion_thruster"].append(commands["ion_thruster"])
            history["xenon"].append(y[IDX_MXE])
            history["wheel_rpm"].append(y[IDX_WHEELS] / wheel_inertia * 60.0 / (2.0 * np.pi))
            history["panel_integrity"].append(float(panel_integrity))
            history["total_impacts"].append(int(total_impacts))
            history["power_proxy_W"].append(float(power_proxy_W))
            history["expected_lambda"].append(float(expected_lambda_tick))
            
            # --- Orbital diagnostics (recording, does not affect physics) ---
            rn_ = np.linalg.norm(y[IDX_POS]); v2_ = np.dot(y[IDX_VEL], y[IDX_VEL])
            a_ = -mu / (2.0 * (0.5*v2_ - mu/rn_))
            h_ = np.cross(y[IDX_POS], y[IDX_VEL])
            i_ = np.degrees(np.arccos(np.clip(h_[2]/np.linalg.norm(h_), -1.0, 1.0)))
            n_vec_ = np.array([-h_[1], h_[0], 0.0])
            Om_ = np.degrees(np.arctan2(n_vec_[1], n_vec_[0])) % 360.0
            e_vec_ = ((v2_ - mu/rn_)*y[IDX_POS] - np.dot(y[IDX_POS], y[IDX_VEL])*y[IDX_VEL]) / mu
            e_ = float(np.linalg.norm(e_vec_))
            eph_d = calculate_body_positions(current_dt, y[IDX_POS])
            sun_ = eph_d["Sun"] / np.linalg.norm(eph_d["Sun"])
            RA_ = np.degrees(np.arctan2(sun_[1], sun_[0])) % 360.0
            dL_ = (Om_ - RA_ + 180.0) % 360.0 - 180.0
            LTAN_ = (12.0 + dL_ / 15.0) % 24.0
            
            history["a_osc"].append(float(a_))
            history["e_osc"].append(e_)
            history["i_osc_deg"].append(float(i_))
            history["LTAN_h"].append(float(LTAN_))

        # 4) Physics: one RK4 step with the current commands
        y, current_dt = rk4_step(y, dt_sim, current_dt, commands)
        t_sim += dt_sim
        tick += 1

        # 5) Pacing WITHOUT drift + listen to ESP32 during the wait (timestamps real message arrivals)
        deadline += cycle_duration
        if deadline - time.perf_counter() <= 0:
            late_cycles += 1
            max_delay = max(max_delay, time.perf_counter() - deadline)
            
        while True:
            wait_time = deadline - time.perf_counter()
            if wait_time <= 0:
                break
            time.sleep(min(wait_time, 0.002))
            msg_queue.extend(connection.read_messages())


    # --- Cleanup live map ---
    if _map_proc is not None:
        try:
            _map_stop.set()
            try:
                _map_queue.put_nowait(None)
            except Exception:
                pass
            _map_proc.join(timeout=2.0)
            if _map_proc.is_alive():
                _map_proc.terminate()
        except Exception:
            pass

    real_duration = time.perf_counter() - t_start
    latencies_ms = np.array(latencies_ms)
    counters = {"bytes_sent": connection.bytes_sent, "bytes_received": connection.bytes_received,
                "skipped_lines": connection.skipped_lines, "ignored_lines": connection.ignored_lines,
                "write_timeouts": connection.write_timeouts, "discarded_fields": discarded_fields,
                "messages_from_esp32": n_msg_esp, "late_cycles": late_cycles}
    
    results = {"history": history, "t_sim": t_sim, "real_duration": real_duration, "effective_k": t_sim / real_duration,
               "late_cycles": late_cycles, "max_delay_ms": max_delay * 1e3, "messages_from_esp32": n_msg_esp,
               "bytes_sent": connection.bytes_sent, "command_delays_s": np.array(cmd_delays),
               "latencies_ms": latencies_ms, "counters": counters,
               "pointing_errors": np.array(pointing_errors), "events": events, "final_state": y, "stop_reason": reason}
               
    if verbose:
        # Print a nice summary at the end of the run. It's always good to have a quick recap.
        print("-" * 60)
        print(f"Finished: {reason}. Simulated {t_sim:.0f} s ({t_sim / 3600:.2f} h) in {real_duration:.1f} real s "
              f"-> effective k {results['effective_k']:.1f} (requested {k:g})")
        print(f"Late cycles: {late_cycles} (max delay {results['max_delay_ms']:.1f} ms)")
        print(f"Messages from ESP32: {n_msg_esp} over {tick} cycles | bytes sent {counters['bytes_sent']} | "
              f"received {counters['bytes_received']}")
        print(f"Skipped lines {counters['skipped_lines']} | discarded fields {discarded_fields} | "
              f"write timeouts {counters['write_timeouts']}")
        if len(latencies_ms):
            print(f"Real round-trip latency: median {np.median(latencies_ms):.1f} ms, "
                  f"p95 {np.percentile(latencies_ms, 95):.1f} ms, max {np.max(latencies_ms):.1f} ms")
        if len(pointing_errors):
            print(f"Camera pointing: max error {np.max(pointing_errors):.4f} deg | mean {np.mean(pointing_errors):.4f} deg")
        print(f"Propellant remaining: xenon {y[IDX_MXE]:.4f} kg | chemical {y[IDX_MCH]:.4f} kg")
        
    return results


def save_results(results, name="hwil_results"):
    # Save the history arrays. Try scipy's .mat format (useful for MATLAB folks),
    # otherwise fallback to standard numpy .npz.
    data = {k: np.array(v) for k, v in results["history"].items()}
    data.update({"latencies_ms": results["latencies_ms"], "command_delays_s": results["command_delays_s"],
                 "pointing_errors": results["pointing_errors"]})
    try:
        from scipy.io import savemat
        savemat(name + ".mat", data)
    except ImportError:
        np.savez(name + ".npz", **data)
        
    with open(name + "_events.json", "w", encoding="utf-8") as f:
        json.dump({"events": results["events"], "counters": results["counters"], "stop_reason": results["stop_reason"]},
                  f, indent=1)


def export_firmware_header(filepath="sim_params.h"):
    """Generates a C header with the constants that firmware and simulator must share."""
    def write_matrix(name, M):
        body = ",\n    ".join("{" + ", ".join(f"{x:.9g}" for x in row) + "}" for row in M)
        return f"static const double {name}[{M.shape[0]}][{M.shape[1]}] = {{\n    {body}\n}};\n"

    lines = ["// AUTOMATICALLY GENERATED by the simulator: DO NOT edit by hand\n", "#pragma once\n\n"]
    
    # Calculate the Julian Date of the initial epoch (needed by the firmware to calculate LTAN)
    from datetime import datetime as _dt
    _dt0 = _dt(year, month, day, hour, minute, second)
    _a = _dt0.year - 1 if _dt0.month <= 2 else _dt0.year
    _m = _dt0.month + 12 if _dt0.month <= 2 else _dt0.month
    _A = int(_a / 100); _B = 2 - _A + int(_A / 4)
    _jd = (int(365.25 * (_a + 4716)) + int(30.6001 * (_m + 1))
           + _dt0.day + _B - 1524.5)
    _jd += (_dt0.hour + _dt0.minute / 60.0 + _dt0.second / 3600.0) / 24.0

    constants = {
        "DT_SIM_S": dt_sim, "WHEEL_LOOP_TAU_S": wheel_loop_tau,
        "JD_EPOCH": _jd,
        "MIN_HALL_FRACTION": min_hall_fraction, "MAX_HALL_THRUST_N": max_hall_thrust,
        "MAX_WHEEL_SPEED_RPM": max_wheel_speed * 60.0 / (2.0 * np.pi),
        "MAX_WHEEL_MOMENTUM_NMS": max_wheel_momentum, "MAX_WHEEL_TORQUE_NM": max_wheel_torque,
        "WHEEL_INERTIA_KGM2": wheel_inertia, "MAX_MTQ_DIPOLE_AM2": max_mtq_dipole,
        "EARTH_MU_KM3S2": mu, "EARTH_EQUATORIAL_RADIUS_KM": earth_equatorial_radius, "EARTH_J2": J2,
    }
    
    for name, value in constants.items():
        lines.append(f"#define {name} {value:.9g}\n")
        
    for name, f in SENSOR_FREQUENCIES_HZ.items():
        lines.append(f"#define FREQ_{name.upper()}_HZ {f:.9g}\n")
        
    for name, comb in RCS_COMBINATIONS.items():
        mask = sum(1 << i for i, b in enumerate(comb) if b)
        lines.append(f"#define RCS_{name.upper().replace('+', 'P').replace('-', 'M')} 0x{mask:02X}\n")
        
    lines.append("\n")
    lines.append(write_matrix("WHEEL_AXES", WHEEL_AXES))
    lines.append(write_matrix("WHEEL_AXES_PINV", WHEEL_AXES_PINV))
    lines.append("static const double WHEEL_NULLSPACE[4] = {0.5, -0.5, 0.5, -0.5};\n")
    lines.append("static const double BODY_INERTIA_KGM2[3] = {"
                 + ", ".join(f"{x:.9g}" for x in np.diag(calculate_inertia_tensor())) + "};\n")
                 
    with open(filepath, "w", encoding="utf-8") as f:
        f.writelines(lines)


# Shared C header: generated directly into the Arduino sketch folder,
# so the firmware finds it right next to "esp32_satellite.ino" (no manual copying needed).
OUT_FIRMWARE = Path.home() / "Desktop" / "hwil" / "firmware" / "esp32_satellite"

# Simulation results: separate folder, outside the Arduino project.
OUT_RESULTS = Path.home() / "Desktop" / "hwil" / "results"

OUT_FIRMWARE.mkdir(parents=True, exist_ok=True)
OUT_RESULTS.mkdir(parents=True, exist_ok=True)

if __name__ == "__main__":
    # Setup a clean Ctrl+C handler.
    # Note: Jupyter often runs on a worker thread and refuses to let us install this.
    # In that case, just use the "Interrupt Kernel" button.
    try:
        signal.signal(signal.SIGINT, _sigint_handler)
    except (ValueError, RuntimeError, AttributeError):
        print("Note: SIGINT handler not installed (Jupyter?). Stop via max_real_time or 'Interrupt'.")

    export_firmware_header(str(OUT_FIRMWARE / "sim_params.h"))
    print(f"Firmware header written to: {OUT_FIRMWARE / 'sim_params.h'}")
    
    connection = SerialConnection()
    results = None
    try:
        results = run_simulation(connection)
    except KeyboardInterrupt:
        print("\nSimulation interrupted by user (Ctrl+C / Interrupt kernel).")
    finally:
        connection.close()
        if results is not None:
            save_results(results, str(OUT_RESULTS / "hwil_results"))
            print(f"Results written to:        {OUT_RESULTS}")

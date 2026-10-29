"""Separate process for rendering the map (cartopy + Tk).

This module is imported by simulator.py and launched via multiprocessing.Process.
It receives a tuple through a queue: (t_sim, lon_sat, lat_sat, sub_lat, sub_lon):
- t_sim: simulated time in seconds, used to trim the ground track to the last 2 orbits.
- sub_lat/sub_lon: subsolar point (where the Sun is at zenith), for the day/night overlay.
"""
import matplotlib

# We force TkAgg here because this runs in a separate process.
# If we let it use the Jupyter inline backend, the window simply won't open.
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np


def map_process(queue, max_points, stop_event, orbit_period_s=5670.0):
    """Target function for the child process. Opens a native window and updates from queue."""
    fig = plt.figure(figsize=(11, 5.5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global()

    # --- Background: land / ocean / coastlines ---
    ax.add_feature(cfeature.LAND, facecolor='#e8e8e8', zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor='#a8c8e8', zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='#555555', zorder=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.2, linestyle=':', alpha=0.4, zorder=1)
    ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.4,
                 xlocs=range(-180, 181, 60), ylocs=range(-90, 91, 30))

    # --- Equator line ---
    ax.plot([-180, 180], [0, 0], color='red', linewidth=1.0,
            linestyle='--', zorder=2, label='Equator')

    # --- Day/Night overlay (RGBA image, updated periodically to save CPU) ---
    lon_grid = np.linspace(-180, 180, 720)
    lat_grid = np.linspace(-90, 90, 360)
    LON, LAT = np.meshgrid(lon_grid, lat_grid)
    
    # Start fully transparent
    rgba = np.zeros((*LON.shape, 4))
    daynight_img = ax.imshow(rgba, extent=[-180, 180, -90, 90],
                             transform=ccrs.PlateCarree(),
                             zorder=0.5, aspect='auto', interpolation='bilinear',
                             origin='lower')

    # --- Ground track and satellite marker ---
    track_line, = ax.plot([], [], color='#ff8800', linewidth=1.2,
                          zorder=3, label='Ground track')
    marker, = ax.plot([], [], 'o', color='red', markersize=9,
                      markeredgecolor='black', markeredgewidth=0.5,
                      zorder=4, label='Satellite')

    ax.legend(loc='upper right', fontsize=8, framealpha=0.8)
    ax.set_title('Sun-Synchronous Orbit — Ground Track', fontsize=10)
    fig.tight_layout()

    plt.ion()
    try:
        plt.show(block=False)
        fig.canvas.draw()
        fig.canvas.flush_events()
    except Exception:
        pass

    # --- Timestamped buffers: ground track limited to the last 2 orbits ---
    t_buf, lon_buf, lat_buf = [], [], []
    last_terminator_update_t = -1e18
    TERMINATOR_UPDATE_INTERVAL_S = 120.0        # [sim s] Update the overlay every 2 sim minutes

    while not stop_event.is_set():
        try:
            # Short timeout to keep the GUI responsive even if no data arrives
            item = queue.get(timeout=0.1)
        except Exception:
            if not plt.fignum_exists(fig.number):
                break
            try:
                fig.canvas.flush_events()
            except Exception:
                pass
            continue

        if item is None:                # Shutdown sentinel
            break

        # Expected format: (t_sim, lon_sat, lat_sat, sub_lat, sub_lon)
        t_sim, lon, lat, sub_lat, sub_lon = item

        # --- Ground track: append new data and drop points older than 2 orbits ---
        t_buf.append(t_sim); lon_buf.append(lon); lat_buf.append(lat)
        t_min = t_sim - 2.0 * orbit_period_s
        while t_buf and t_buf[0] < t_min:
            t_buf.pop(0); lon_buf.pop(0); lat_buf.pop(0)
        while len(t_buf) > max_points:   # Safety net just in case
            t_buf.pop(0); lon_buf.pop(0); lat_buf.pop(0)

        # Split the ground track at the antimeridian (|Δlon| > 180).
        # Classic matplotlib hack: inserting NaNs breaks the line so it doesn't wrap around the globe.
        if len(lon_buf) < 2:
            lons, lats = list(lon_buf), list(lat_buf)
        else:
            lons = [lon_buf[0]]; lats = [lat_buf[0]]
            for i in range(1, len(lon_buf)):
                if abs(lon_buf[i] - lon_buf[i-1]) > 180.0:
                    lons.append(float('nan')); lats.append(float('nan'))
                lons.append(lon_buf[i]); lats.append(lat_buf[i])

        marker.set_data([lon], [lat])
        track_line.set_data(lons, lats)
        ax.set_title(
            f'Sun-Synchronous Orbit — LAT {lat:+06.2f}°  LON {lon:+07.2f}°  '
            f'(t={t_sim/3600:.2f} h)', fontsize=10)

        # --- Day/Night overlay: update every TERMINATOR_UPDATE_INTERVAL_S sim seconds ---
        if t_sim - last_terminator_update_t >= TERMINATOR_UPDATE_INTERVAL_S:
            last_terminator_update_t = t_sim
            
            # Compute solar zenith angle using spherical trigonometry.
            # cos(zenith) = sin(lat)sin(sub_lat) + cos(lat)cos(sub_lat)cos(lon-sub_lon)
            lat_r = np.radians(LAT)
            sub_lat_r = np.radians(sub_lat)
            cos_z = (np.sin(lat_r) * np.sin(sub_lat_r) +
                     np.cos(lat_r) * np.cos(sub_lat_r) *
                     np.cos(np.radians(LON - sub_lon)))
            
            # Night (cos_z < 0): dark alpha; Day: transparent.
            # We apply a smoothstep on a ±BAND_WIDTH transition to avoid jagged edges.
            BAND_WIDTH = 0.05                              # ≈ 3° of penumbra on each side
            u = np.clip((cos_z + BAND_WIDTH) / (2.0 * BAND_WIDTH), 0.0, 1.0)   # 0=night, 1=day
            smooth = u * u * (3.0 - 2.0 * u)               # Smoothstep polynomial
            rgba[..., :3] = 0.03                           # Almost black
            rgba[..., 3] = (1.0 - smooth) * 0.45           # Night=0.45 alpha, Day=0
            daynight_img.set_data(rgba)

        if not plt.fignum_exists(fig.number):
            break

        try:
            fig.canvas.draw_idle()
            fig.canvas.flush_events()
        except Exception:
            pass

    try:
        plt.close('all')
    except Exception:
        pass

# HWIL Satellite Simulator: Sun-Synchronous LEO Mission

![Python](https://img.shields.io/badge/Python-3.10-blue)
![C++](https://img.shields.io/badge/C++-Arduino-blue)
![Platform](https://img.shields.io/badge/Platform-ESP32-orange)
![License](https://img.shields.io/badge/License-MIT-green)

## Overview

This project is a **Hardware-In-The-Loop (HWIL) Digital Twin** for a 150 kg microsatellite in a Sun-Synchronous Orbit (SSO) at ~500 km altitude (Sentinel-2 style).

The goal is to test the real flight software—running on an **ESP32** microcontroller—against a high-fidelity physics simulation running on a host PC. The PC simulates the space environment and the sensors, while the ESP32 acts as the on-board computer (OBC), making autonomous decisions based on the incoming telemetry.

## Key Features

*   **High-Fidelity Orbital Mechanics**: RK4 integrator with J2 perturbation, atmospheric drag (NRLMSISE-00 simplified), Solar Radiation Pressure (SRP), and Third-Body perturbations (Sun & Moon).
*   **Attitude Dynamics (ADCS)**: Quaternion-based kinematics and Euler dynamics. PD controller acting on a 4-wheel pyramidal reaction wheel assembly.
*   **Autonomous Orbit Maintenance**: The ESP32 maintains the SSO and Repeat Ground Track (RGT) by calculating required \(\Delta v\) and firing a Hall-effect thruster (15 mN class).
*   **Space Debris & Collision Avoidance**:
    *   *Small debris*: Poisson-distributed impacts that degrade solar panel efficiency over time.
    *   *Large debris*: Real-time Conjunction Data Messages (CDMs). The ESP32 autonomously plans and executes Collision Avoidance Maneuvers (CAM) when a warning is received.
*   **Embedded I/O & Hardware Validation**:
    *   **I/O Loopback Self-Test**: A hardware verification of the digital I/O channel (GPIO 18 → GPIO 19) to ensure interface integrity during the entire run.
    *   **Status LED**: A physical LED (GPIO 2) provides real-time visual feedback of the Hall thruster state, confirming the end-to-end chain from simulation to actuator command.
*   **Live Ground Track**: A multiprocessing 2D map (Cartopy) showing the satellite's real-time position, ground track history, and the day/night terminator.
*   **ASCII Dashboard**: A real-time terminal dashboard displaying orbital elements, propellant levels, actuator states, and serial latency metrics.

## Testing & Validation

This project includes a formal test campaign to validate the HWIL setup, including interface integrity, I/O loopback, load testing, and end-to-end chain validation.

📄 **The complete test procedures and results are documented in [TEST_PROCEDURES.md](TEST_PROCEDURES.md).**

## System Architecture

The system is split into two halves that communicate over a serial link using line-delimited JSON at 115200 baud.

```text
+-----------------------------+         Serial / JSON         +-----------------------------+
|        HOST PC (Python)     | <---------------------------> |      OBC (ESP32 / C++)      |
|                             |                               |                             |
|  - Physics Engine (RK4)     |   PC sends: Sensor Data       |   - Attitude Control (PD)   |
|  - Sensor Emulation         |   ESP32 sends: Actuator Cmds  |   - Orbit Maintenance       |
|  - Space Debris Models      |                               |   - CAM Planner             |
|  - 2D Map Rendering         |                               |   - I/O Self-Test & LED     |
+-----------------------------+                               +-----------------------------+
```

## Repository Structure

```text
hwil/
├── simulator/                 # Python simulation environment
│   ├── simulator.py           # Core physics, sensors, and communication loop
│   ├── map_process.py         # Multiprocessing script for the live map
│   ├── main.py                # Launcher script
│   └── requirements.txt       # Python dependencies
├── firmware/
│   └── esp32_satellite/
│       └── esp32_satellite.ino # Flight software for the ESP32 (Arduino IDE)
├── README.md
├── TEST_PROCEDURES.md         # Formal test procedures and results
└── .gitignore
```

## Getting Started

### Prerequisites

*   **Python 3.10+**
*   **Arduino IDE** (with the ESP32 board manager installed)
*   **ArduinoJson** library (v7.x) installed via Arduino IDE Library Manager.
*   **Hardware**: An ESP32 development board (e.g., ESP32 DevKitC) connected via USB.
*   **Hardware Setup**: A jumper wire connecting **GPIO 18** to **GPIO 19** (for the I/O loopback self-test).

### Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/A-Sasha-A/satellite-hwil-simulator.git
    cd satellite-hwil-simulator
    ```

2.  Install the Python dependencies:
    ```bash
    pip install -r simulator/requirements.txt
    ```

### Hardware Setup & Firmware Flash

1.  Connect a jumper wire between **GPIO 18** and **GPIO 19** on the ESP32.
2.  Open `firmware/esp32_satellite/esp32_satellite.ino` in the Arduino IDE.
3.  The firmware requires the `sim_params.h` file. Generate it by running the Python script once (see step 1 of "Running the Simulation" below).
4.  Select the correct board (`ESP32 Dev Module`) and port in Arduino IDE.
5.  Compile and upload the firmware to the ESP32.
6.  **Important**: After the upload, **close the Arduino Serial Monitor**. The port must be free for Python to use.

## Running the Simulation

You can launch the simulation from the Terminal or via Jupyter Notebook.

### Option 1: Terminal (Recommended)

1.  Open a terminal and navigate to the simulator folder:
    ```bash
    cd simulator
    ```
2.  Run the launcher:
    ```bash
    python main.py
    ```

### Option 2: Jupyter Notebook

Create a notebook inside the `simulator/` folder and run the following cell:

```python
import sys, os
from pathlib import Path
os.chdir(Path.home() / "Desktop" / "hwil" / "simulator")
sys.path.insert(0, str(Path.cwd()))

import importlib
import simulator
importlib.reload(simulator)

# Generate the firmware header (if not already done)
simulator.export_firmware_header(str(simulator.OUT_FIRMWARE / "sim_params.h"))

# Start the simulation
link = simulator.SerialConnection()
try:
    results = simulator.run_simulation(link)
finally:
    link.close()
    if results is not None:
        simulator.save_results(results, str(simulator.OUT_RESULTS / "hwil_results"))
```

## Physics & Control Models

### Environmental Perturbations
*   **Orbit**: The state vector is propagated using a fixed-step Runge-Kutta 4 (RK4) integrator.
*   **J2**: Accounts for Earth's oblateness, critical for maintaining the SSO node drift.
*   **Atmosphere**: Exponential density model scaled by a solar activity factor (F10.7).
*   **SRP**: Solar Radiation Pressure calculated by projecting the Sun vector onto the satellite's flat surfaces (optical properties included).

### Attitude Determination & Control
The ESP32 receives quaternion data from a virtual star tracker and angular velocity from a virtual gyroscope. It runs a PD controller to keep the camera (body +z) pointing at the nadir. Actuation is provided by four reaction wheels in a pyramidal configuration. Excess momentum is dumped using three orthogonal magnetorquers (\(\tau = m \times B\)).

### Collision Avoidance (CAM)
When the Python script detects a high-probability conjunction, it sends a warning (`conj`) with the Time of Closest Approach (TCA) and a target geocentric radius. The ESP32 evaluates the risk and, if time permits, fires the Hall thruster along-track (body +x) to raise the orbit and increase the miss distance. If the CAM is successful, the orbital radius will be greater than or equal to the target at the TCA.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Alessandro Antonio Micali Capparelli**
Aerospace Engineering Student / Researcher
*Developed as a practical demonstration of HWIL simulation for autonomous spacecraft operations.*

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

import sys
from pathlib import Path

# Make sure the custom simulator folder is in the system path.
# (Hardcoding paths is not ideal for production, but it's fine for this local HWIL setup).
sys.path.insert(0, str(Path.home() / "Desktop" / "hwil" / "simulator"))

# Import the core simulator.
# Note: do NOT paste the simulator code directly into a Jupyter cell,
# multiprocessing on macOS will throw a fit and the map process won't start.
import simulator
import importlib

# Reload just in case we are testing inside Jupyter and modified the source code
importlib.reload(simulator)

def main():
    # Open the serial connection to the ESP32 board
    link = simulator.SerialConnection()
    results = None
    
    try:
        # Start the main loop and let the physics do its thing
        results = simulator.run_simulation(link)
    except KeyboardInterrupt:
        print("\n[!] Simulation interrupted by user. Closing down cleanly...")
    finally:
        # Always close the serial port, even if the simulation crashes half-way through.
        # Otherwise the port stays locked and you have to unplug the USB cable.
        link.close()
        
        # If we actually got some data, save the telemetry and event logs to disk
        if results is not None:
            simulator.save_results(results, str(simulator.OUT_RESULTS / "hwil_results"))
            print(f"Results successfully written to: {simulator.OUT_RESULTS}")

if __name__ == "__main__":
    main()

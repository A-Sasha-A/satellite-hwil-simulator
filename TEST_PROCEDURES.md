# HWIL Test Procedures and Results

## Test Environment
- Host: macOS (Apple Silicon)
- Python: 3.10
- On-Board Computer: ESP32 DevKit (dual-core, 240 MHz)
- Serial link: USB-UART (CH340), 115200 baud, JSON line protocol
- I/O Loopback: GPIO 18 (OUT) → GPIO 19 (IN)

## Test 1 — Serial Interface Integrity
**Objective**: Verify that the PC↔ESP32 serial link is functional.
**Procedure**: Send `{"cmd":"ping"}`, wait for `{"log":"pong"}` within 100 ms, repeat 100 times.
**Pass criteria**: ≥99% pong replies, no corrupted lines.
**Result**: PASSED — 100/100 replies, median RTT 30.4 ms, p95 49.9 ms.

## Test 2 — GPIO I/O Loopback Self-Test
**Objective**: Verify the integrity of the digital I/O channel.
**Procedure**: Drive GPIO 18 with a known pattern (LOW/HIGH), read back on GPIO 19 with a 200 µs propagation delay, compare.
**Pass criteria**: 4/4 patterns match.
**Result**: PASSED — 4/4 after increasing delay from 20 µs to 200 µs (breadboard capacitance issue found and resolved).

## Test 3 — Load Test (Sensor Throughput)
**Objective**: Verify the ESP32 can process sensor packets without buffer overflow.
**Procedure**: Increase GPS sensor frequency from 0.1 Hz to 1 Hz simulated, run for 10 minutes, monitor `Skipped lines` and `Discarded fields` in the dashboard.
**Pass criteria**: Zero skipped lines, zero discarded fields.
**Result**: PASSED — 0 errors over 600+ packets.

## Test 4 — End-to-End Chain Test (Collision Avoidance)
**Objective**: Verify the complete chain: simulation → OBC decision → actuator hardware.
**Procedure**: Inject a conjunction warning, verify the ESP32 commands Hall thruster, verify the LED turns on, verify the CAM completes and the LED turns off.
**Pass criteria**: LED ON during maneuver, LED OFF after completion, no I/O warnings.
**Result**: PASSED — LED responded correctly to Hall state changes.

## Known Issues and Resolutions
- **Issue**: I/O loopback self-test failed intermittently with 20 µs delay.
  **Root Cause**: Breadboard capacitance caused signal propagation delay longer than expected.
  **Resolution**: Increased delay to 200 µs. No further failures.

## Conclusion
The HWIL setup successfully validates the flight software's behavior against a high-fidelity spacecraft simulator. All critical interfaces have been tested and verified.

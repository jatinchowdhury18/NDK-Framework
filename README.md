# Nodal-DK Code Generator

This repository contains a system for generating C++/Eigen
code for computing circuit simulations using the nodal DK
method. The workflow is as follows:
- Write a JSON circuit configuration (see provided examples)
- Run the `generate_ndk_cpp.py` script to generate the C++ code
- Compile and run your C++ code

Disclaimer: this system has not been rigorously tested,
use at your own risk! If you run into any issues feel free
to make an issue, but know that I may or may not have time
to actually fix the issue.

## Notes:

Undocumented netlist gotchas:
- Ground is node 0
- Netlist entries go from positive node to negative
- Diodes go from anode to cathode
- Transistors go "Base, Collector, Emitter"
- Currently only supports one "input" and one "output"

Supported elements:
- Resistor (R)
- Capacitor (C)
- Inductor (L)
- Voltage Source (V)
- Diode (D)
- Diode Pair (DP)
- Transistor (Q)
- Ideal Op-Amp (A)

Unsupported elements:
- Triode (WIP)
- OTA

Example circuits:
- [RC Lowpass Filter](https://www.electronics-tutorials.ws/filter/filter_2.html)
- [Op-Amp Sallen-Key Lowpass Filter](https://www.electronics-tutorials.ws/filter/sallen-key-filter.html)
- [Cry Baby Circuit](https://www.electrosmash.com/crybaby-gcb-95)
- [Modified Big Muff Pi Clipping Circuit](https://www.electrosmash.com/big-muff-pi-analysis)

## License

This software is provided under th MIT license.

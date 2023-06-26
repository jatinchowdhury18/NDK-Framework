# Nodal-DK Code Generator

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

Unsupported elements:
- Triode
- Op-Amp
- OTA

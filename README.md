# Nodal-DK Code Generator

Undocumented netlist gotchas:
- Ground is node 0
- Netlist entried go from positive node to negative
- Diodes go from anode to cathode
- Transistors go "Base, Collector, Emitter"
- Variable voltage sources must come before fixed voltgae sources
- Currently only supports one "input" and one "output"

Supported elements:
- Resistor
- Capacitor
- Inductor
- Voltage Source
- Diode
- Diode Pair
- Transistor

Unsupported elements:
- Triode
- Op-Amp

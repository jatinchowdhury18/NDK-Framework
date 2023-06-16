from dataclasses import dataclass
from enum import Enum

class Element(Enum):
    Resistor = 1
    Capacitor = 2
    Inductor = 3
    Voltage = 4
    Diode = 5 # TODO!
    Transistor = 6

ONE_PORT_ELEMENTS = [Element.Resistor, Element.Capacitor, Element.Voltage, Element.Inductor]

@dataclass
class ElementInfo:
    nodes: list[int]
    name: str = ""
    element: Element = Element.Resistor
    is_constant: bool = True

class NetlistInfo:
    num_nodes: int = 0
    num_resistors: int = 0
    num_voltages: int = 0
    num_pots: int = 0
    num_states: int = 0
    num_nl_ports: int = 0
    elements: list[ElementInfo] = []

def parse_netlist(netlist: list[str]):
    def get_element_type(tag: str) -> Element:
        if tag[0] == 'R':
            return Element.Resistor
        elif tag[0] == 'C':
            return Element.Capacitor
        elif tag[0] == 'L':
            return Element.Inductor
        elif tag[0] == 'V':
            return Element.Voltage
        elif tag[0] == 'Q':
            return Element.Transistor
        assert False, f"Unknown element type!"

    info = NetlistInfo()
    for line in netlist:
        line_parts = line.split(' ')
        el_type = get_element_type(line_parts[0])
        print(f"Adding element: {el_type}")

        if el_type in ONE_PORT_ELEMENTS:
            el = ElementInfo(nodes=[int(line_parts[1]), int(line_parts[2])])
            el.name = line_parts[0]
            el.element = el_type

            if el_type is Element.Resistor and ('VARIABLE' in line_parts):
                el.is_constant = False
            elif el_type is Element.Voltage and not ('FIXED' in line_parts):
                el.is_constant = False
        elif el_type is Element.Transistor:
            el = ElementInfo(nodes=[int(line_parts[1]), int(line_parts[2]), int(line_parts[3])])
            el.name = line_parts[0]
            el.element = el_type

        info.elements.append(el)

        info.num_nodes = max(info.num_nodes, max(el.nodes))
        if el.element is Element.Resistor:
            if el.is_constant:
                info.num_resistors += 1
            else:
                info.num_pots += 1
        if el.element is Element.Voltage:
            info.num_voltages += 1
        if el.element in [Element.Capacitor, Element.Inductor]:
            info.num_states += 1
        if el.element is Element.Transistor:
            info.num_nl_ports += 2
        

    print(f"Circuit contains: {info.num_nodes} nodes, {info.num_voltages} voltage sources, {info.num_resistors} resistors, {info.num_states} states, {info.num_pots} pots, {info.num_nl_ports} nonlinear ports")

    return info

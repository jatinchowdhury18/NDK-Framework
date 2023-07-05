from dataclasses import dataclass
from enum import Enum

class Element(Enum):
    Resistor = 1
    Capacitor = 2
    Inductor = 3
    Voltage = 4
    Diode = 5
    DiodePair = 6
    Transistor = 7

ONE_PORT_ELEMENTS = [Element.Resistor,
                     Element.Capacitor,
                     Element.Voltage,
                     Element.Inductor,
                     Element.Diode,
                     Element.DiodePair]

@dataclass
class ElementInfo:
    nodes: list[int]
    name: str = ""
    element: Element = Element.Resistor
    is_constant: bool = True

class NetlistInfo:
    num_nodes: int = 0
    num_resistors: int = 0
    num_var_voltages: int = 0
    num_fixed_voltages: int = 0
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
        elif tag[0:2] == 'DP':
            return Element.DiodePair
        elif tag[0] == 'D':
            return Element.Diode
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

        if el.element is Element.Voltage and not el.is_constant:
            insert_idx = 0
            for idx, el_test in enumerate(info.elements):
                if el_test is not Element.Voltage or el_test.is_constant:
                    insert_idx = idx + 1
                    break
            info.elements.insert(insert_idx, el)
        else:
            info.elements.append(el)

        info.num_nodes = max(info.num_nodes, max(el.nodes))
        if el.element is Element.Resistor:
            if el.is_constant:
                info.num_resistors += 1
            else:
                info.num_pots += 1
        if el.element is Element.Voltage:
            if el.is_constant:
                info.num_fixed_voltages += 1
            else:
                info.num_var_voltages += 1
            info.num_voltages += 1
        if el.element in [Element.Capacitor, Element.Inductor]:
            info.num_states += 1
        if el.element is Element.Transistor:
            info.num_nl_ports += 2
        if el.element is Element.Diode or el.element is Element.DiodePair:
            info.num_nl_ports += 1
        

    print(f"Circuit contains: {info.num_nodes} nodes, {info.num_voltages} voltage sources, {info.num_resistors} resistors, {info.num_states} states, {info.num_pots} pots, {info.num_nl_ports} nonlinear ports")

    return info

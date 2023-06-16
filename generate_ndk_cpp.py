import json
import sys
from utils.netlist_parser import parse_netlist
from utils.generate_header import generate_header
from utils.generate_cpp import generate_cpp

config_file = sys.argv[1]
print(f"Generating NDK simulator for: {config_file}")

with open(config_file) as f:
    config_json = json.load(f)

netlist_info = parse_netlist(config_json['netlist'])
# print(netlist_info.elements)
circuit_outputs = config_json['output_nodes']
num_outputs = len(circuit_outputs)


struct_name = config_json['struct_name']
header_file = generate_header(config_json, netlist_info, num_outputs)
with open(f'{struct_name}.h', 'w') as f:
    for line in header_file:
        f.write(line + '\n')

cpp_file = generate_cpp(config_json, netlist_info, circuit_outputs)
with open(f'{struct_name}.cpp', 'w') as f:
    for line in cpp_file:
        f.write(line + '\n')

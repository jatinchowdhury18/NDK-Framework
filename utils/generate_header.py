import sys
from datetime import datetime
from .netlist_parser import NetlistInfo

def generate_header(config_json, netlist_info: NetlistInfo, num_outputs):
    header_file = []

    header_file.append("/*")
    header_file.append(f" * This file was generated on {datetime.now()}")
    header_file.append(f" * using the command: `{' '.join(sys.argv)}`")
    header_file.append(" */")

    header_file.append("#pragma once\n")
    header_file.append("// START USER INCLUDES")
    for include_line in config_json['header_includes']:
        header_file.append(include_line)
    header_file.append("// END USER INCLUDES\n")

    if 'cpp_namespace' in config_json:
        header_file.append(f"namespace {config_json['cpp_namespace']} {{")

    struct_name = config_json['struct_name']
    header_file.append(f"struct {struct_name}")
    header_file.append("{")

    header_file.append("    // START USER ENTRIES")
    for include_line in config_json['cpp_struct_entries']:
        header_file.append(include_line)
    header_file.append("    // END USER ENTRIES\n")

    header_file.append(f"    using T = double;")
    header_file.append(f"    static constexpr int num_nodes = {netlist_info.num_nodes};")
    header_file.append(f"    static constexpr int num_resistors = {netlist_info.num_resistors};")
    header_file.append(f"    static constexpr int num_pots = {netlist_info.num_pots};")
    header_file.append(f"    static constexpr int num_states = {netlist_info.num_states};")
    header_file.append(f"    static constexpr int num_nl_ports = {netlist_info.num_nl_ports};")
    header_file.append(f"    static constexpr int num_voltages = {netlist_info.num_voltages};")
    header_file.append(f"    static constexpr int num_outputs = {num_outputs};")
    header_file.append(f"")

    header_file.append(f"    // State Variables")
    header_file.append(f"    std::array<Eigen::Vector<T, num_states>, MAX_NUM_CHANNELS> x_n;")
    header_file.append(f"    std::array<Eigen::Vector<T, num_nl_ports>, MAX_NUM_CHANNELS> v_n;")
    header_file.append("")

    header_file.append(f"    // NDK Matrices")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_states> A_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_voltages> B_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_nl_ports> C_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_states> D_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_voltages> E_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_nl_ports> F_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_states> G_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_voltages> H_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_nl_ports> K_mat;")
    header_file.append("")

    header_file.append(f"    // Intermediate matrices for NDK matrix updates;")
    header_file.append(f"    Eigen::Matrix<T, num_pots, num_pots> Q;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_pots> Ux;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_pots> Uo;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_pots> Un;")
    header_file.append(f"    Eigen::Matrix<T, num_voltages, num_pots> Uu;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_states> A0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_voltages> B0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_nl_ports> C0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_states> D0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_voltages> E0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_outputs, num_nl_ports> F0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_states> G0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_voltages> H0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_nl_ports, num_nl_ports> K0_mat;")
    header_file.append(f"    Eigen::Matrix<T, num_states, num_states> two_Z_Gx;")
    header_file.append("")

    header_file.append(f"    void reset (T fs);")
    header_file.append(f"    void update_pots (const std::array<T, num_pots>& pot_values);")
    process_dtype = config_json['process_data_type'] if 'process_data_type' in config_json else "T"
    header_file.append(f"    void process (std::span<{process_dtype}> channel_data, size_t channel_index) noexcept;")

    header_file.append("};")

    if 'cpp_namespace' in config_json:
        header_file.append(f"}} // namespace {config_json['cpp_namespace']}")

    return header_file

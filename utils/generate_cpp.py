import sys
from datetime import datetime
from .netlist_parser import *

def generate_Gr(netlist_info):
    Gr_entries = []
    for el in netlist_info.elements:
        if el.element is Element.Resistor and el.is_constant:
            Gr_entries.append(f"(T) 1 / {el.name}")
    return f"    const auto Gr = Eigen::DiagonalMatrix<T, num_resistors> {{ {', '.join(Gr_entries)} }};"

def generate_Gx_Z(netlist_info):
    Gx_entries = []
    Z_entries = []
    for el in netlist_info.elements:
        if el.element is Element.Capacitor:
            Gx_entries.append(f"(T) 2 * fs * {el.name}")
            Z_entries.append("1")
        elif el.element is Element.Inductor:
            Gx_entries.append(f"(T) 1 / ((T) 2 * fs * {el.name})")
            Z_entries.append("-1")
    Gx =  f"    const auto Gx = Eigen::DiagonalMatrix<T, num_states> {{ {', '.join(Gx_entries)} }};"
    Z =  f"    const auto Z = Eigen::DiagonalMatrix<T, num_states> {{ {', '.join(Z_entries)} }};"
    return f"{Gx}\n{Z}\n"

def generate_N_matrices(netlist_info: NetlistInfo, outputs):
    lines = []
    lines.append(f"    // Set up component-defined matrices")

    # Nr matrix
    lines.append(f"    const Eigen::Matrix<T, num_resistors, num_nodes> Nr {{")
    for el in netlist_info.elements:
        if el.element is Element.Resistor and el.is_constant:
            node_map = ["+0"] * netlist_info.num_nodes
            if el.nodes[0] > 0: node_map[el.nodes[0] - 1] = "+1" # positive node
            if el.nodes[1] > 0: node_map[el.nodes[1] - 1] = "-1" # negative node
            lines.append(f"        {{ {', '.join(node_map)} }},")
    lines.append("    };\n")

    # Nv matrix
    if netlist_info.num_pots > 0:
        lines.append(f"    const Eigen::Matrix<T, num_pots, num_nodes> Nv {{")
        for el in netlist_info.elements:
            if el.element is Element.Resistor and not el.is_constant:
                node_map = ["+0"] * netlist_info.num_nodes
                if el.nodes[0] > 0: node_map[el.nodes[0] - 1] = "+1" # positive node
                if el.nodes[1] > 0: node_map[el.nodes[1] - 1] = "-1" # negative node
                lines.append(f"        {{ {', '.join(node_map)} }},")
        lines.append("    };\n")

    # Nx matrix
    lines.append(f"    const Eigen::Matrix<T, num_states, num_nodes> Nx {{")
    for el in netlist_info.elements:
        if el.element in [Element.Capacitor, Element.Inductor]:
            node_map = ["+0"] * netlist_info.num_nodes
            if el.nodes[0] > 0: node_map[el.nodes[0] - 1] = "+1" # positive node
            if el.nodes[1] > 0: node_map[el.nodes[1] - 1] = "-1" # negative node
            lines.append(f"        {{ {', '.join(node_map)} }},")
    lines.append("    };\n")

    # Nu matrix
    lines.append(f"    const Eigen::Matrix<T, num_voltages, num_nodes> Nu {{")
    for el in netlist_info.elements:
        if el.element is Element.Voltage:
            node_map = ["+0"] * netlist_info.num_nodes
            if el.nodes[0] > 0: node_map[el.nodes[0] - 1] = "+1" # positive node
            if el.nodes[1] > 0: node_map[el.nodes[1] - 1] = "-1" # negative node
            lines.append(f"        {{ {', '.join(node_map)} }},")
    lines.append("    };\n")

    # Nn matrix
    lines.append(f"    const Eigen::Matrix<T, num_nl_ports, num_nodes> Nn {{")
    for el in netlist_info.elements:
        if el.element is Element.Transistor:
            base_node_idx = el.nodes[0] - 1
            collector_node_idx = el.nodes[1] - 1
            emitter_node_idx = el.nodes[2] - 1

            node_map_cb = ["+0"] * netlist_info.num_nodes
            if el.nodes[1] > 0: node_map_cb[collector_node_idx] = "+1" # positive node
            if el.nodes[0] > 0: node_map_cb[base_node_idx] = "-1" # negative node
            lines.append(f"        {{ {', '.join(node_map_cb)} }},")

            node_map_ce = ["+0"] * netlist_info.num_nodes
            if el.nodes[1] > 0: node_map_ce[collector_node_idx] = "+1" # positive node
            if el.nodes[2] > 0: node_map_ce[emitter_node_idx] = "-1" # negative node
            lines.append(f"        {{ {', '.join(node_map_ce)} }},")
        elif el.element is Element.Diode or el.element is Element.DiodePair:
            node_map = ["+0"] * netlist_info.num_nodes
            if el.nodes[0] > 0: node_map[el.nodes[0] - 1] = "+1" # positive node
            if el.nodes[1] > 0: node_map[el.nodes[1] - 1] = "-1" # negative node
            lines.append(f"        {{ {', '.join(node_map)} }},")

    lines.append("    };\n")

    # No matrix 
    num_outputs = len(outputs)
    lines.append(f"    const Eigen::Matrix<T, num_outputs, num_nodes> No {{")
    for node_idx in outputs:
        node_map = ["+0"] * netlist_info.num_nodes
        node_map[node_idx - 1] = "+1"
        lines.append(f"        {{ {', '.join(node_map)} }},")
    lines.append("    };\n")

    lines.append(f"    static constexpr auto S_dim = num_nodes + num_voltages;")

    # Nx_0
    lines.append(f"    Eigen::Matrix<T, num_states, S_dim> Nx_0 = Eigen::Matrix<T, num_states, S_dim>::Zero();",)
    lines.append(f"    Nx_0.block<num_states, num_nodes> (0, 0) = Nx;",)

    # Nn_0
    lines.append(f"    Eigen::Matrix<T, num_nl_ports, S_dim> Nn_0 = Eigen::Matrix<T, num_nl_ports, S_dim>::Zero();",)
    lines.append(f"    Nn_0.block<num_nl_ports, num_nodes> (0, 0) = Nn;",)

    # No_0
    lines.append(f"    Eigen::Matrix<T, num_outputs, S_dim> No_0 = Eigen::Matrix<T, num_outputs, S_dim>::Zero();",)
    lines.append(f"    No_0.block<num_outputs, num_nodes> (0, 0) = No;",)

    # Nv_0
    if netlist_info.num_pots > 0:
        lines.append(f"    Eigen::Matrix<T, num_pots, S_dim> Nv_0 = Eigen::Matrix<T, num_pots, S_dim>::Zero();",)
        lines.append(f"    Nv_0.block<num_pots, num_nodes> (0, 0) = Nv;",)

    # Zero_I
    lines.append(f"    Eigen::Matrix<T, num_voltages, S_dim> Zero_I = Eigen::Matrix<T, num_voltages, S_dim>::Zero();",)
    lines.append(f"    Zero_I.block<{netlist_info.num_voltages}, {netlist_info.num_voltages}> (0, {netlist_info.num_nodes}).setIdentity();")
    lines.append("")

    # S0
    lines.append(f"    Eigen::Matrix<T, S_dim, S_dim> S0_mat = Eigen::Matrix<T, S_dim, S_dim>::Zero();")
    lines.append(f"    S0_mat.block<num_nodes, num_nodes> (0, 0) = Nr.transpose() * Gr * Nr + Nx.transpose() * Gx * Nx;")
    lines.append(f"    S0_mat.block<num_nodes, num_voltages> (0, num_nodes) = Nu.transpose();")
    lines.append(f"    S0_mat.block<num_voltages, num_nodes> (num_nodes, 0) = Nu;")
    lines.append(f"    Eigen::Matrix<T, S_dim, S_dim> S0_inv = S0_mat.inverse();\n")

    return lines

def generate_NDK_matrices(netlist_info: NetlistInfo):
    lines = []
    lines.append(f"    // Pre-compute NDK and intermediate matrices")

    if netlist_info.num_pots > 0:
        lines.append("    Q = Nv_0 * (S0_inv * Nv_0.transpose());")
        lines.append("    Ux = Nx_0 * (S0_inv * Nv_0.transpose());")
        lines.append("    Uo = No_0 * (S0_inv * Nv_0.transpose());")
        lines.append("    Un = Nn_0 * (S0_inv * Nv_0.transpose());")
        lines.append("    Uu = Zero_I * (S0_inv * Nv_0.transpose());")
        lines.append("    A0_mat = (T) 2 * (Z * (Gx * (Nx_0 * (S0_inv * Nx_0.transpose())))) - Z.toDenseMatrix();")
        lines.append("    B0_mat = (T) 2 * (Z * (Gx * (Nx_0 * (S0_inv * Zero_I.transpose()))));")
        lines.append("    C0_mat = (T) 2 * (Z * (Gx * (Nx_0 * (S0_inv * Nn_0.transpose()))));")
        lines.append("    D0_mat = No_0 * (S0_inv * Nx_0.transpose());")
        lines.append("    E0_mat = No_0 * (S0_inv * Zero_I.transpose());")
        lines.append("    F0_mat = No_0 * (S0_inv * Nn_0.transpose());")
        lines.append("    G0_mat = Nn_0 * (S0_inv * Nx_0.transpose());")
        lines.append("    H0_mat = Nn_0 * (S0_inv * Zero_I.transpose());")
        lines.append("    K0_mat = Nn_0 * (S0_inv * Nn_0.transpose());")
        lines.append("    two_Z_Gx = (T) 2 * (Z.toDenseMatrix() * Gx.toDenseMatrix());\n")
    else:
        u_fixed = []
        for el in netlist_info.elements:
            if el.element is Element.Voltage:
                if el.is_constant:
                    u_fixed.append(el.name)
        u_fixed_str = ', '.join(u_fixed)

        lines.append("    A_mat = (T) 2 * (Z * (Gx * (Nx_0 * (S0_inv * Nx_0.transpose())))) - Z.toDenseMatrix();")
        lines.append("    Eigen::Matrix<T, num_states, num_voltages> B_mat = (T) 2 * (Z * (Gx * (Nx_0 * (S0_inv * Zero_I.transpose()))));")
        lines.append("    B_mat_var = B_mat.leftCols<num_voltages_variable>();")
        lines.append(f"    B_u_fix = B_mat.rightCols<num_voltages_constant>() * Eigen::Vector<T, num_voltages_constant> {{ {u_fixed_str} }};")
        lines.append("    C_mat = (T) 2 * (Z * (Gx * (Nx_0 * (S0_inv * Nn_0.transpose()))));")
        lines.append("    D_mat = No_0 * (S0_inv * Nx_0.transpose());")
        lines.append("    Eigen::Matrix<T, num_outputs, num_voltages> E_mat = No_0 * (S0_inv * Zero_I.transpose());")
        lines.append("    E_mat_var = E_mat.leftCols<num_voltages_variable>();")
        lines.append(f"    E_u_fix = E_mat.rightCols<num_voltages_constant>() * Eigen::Vector<T, num_voltages_constant> {{ {u_fixed_str} }};")
        lines.append("    F_mat = No_0 * (S0_inv * Nn_0.transpose());")
        lines.append("    G_mat = Nn_0 * (S0_inv * Nx_0.transpose());")
        lines.append("    Eigen::Matrix<T, num_nl_ports, num_voltages> H_mat = Nn_0 * (S0_inv * Zero_I.transpose());")
        lines.append("    H_mat_var = H_mat.leftCols<num_voltages_variable>();")
        lines.append(f"    H_u_fix = H_mat.rightCols<num_voltages_constant>() * Eigen::Vector<T, num_voltages_constant> {{ {u_fixed_str} }};")
        lines.append("    K_mat = Nn_0 * (S0_inv * Nn_0.transpose());\n")


    return lines

def generate_reset_states(config_json):
    lines = []
    lines.append("    // reset state vectors")
    lines.append("    for (size_t ch = 0; ch < MAX_NUM_CHANNELS; ++ch)\n    {")
    if 'initial_state_x_n' in config_json:
        lines.append(f"        x_n[ch] = Eigen::Vector<T, num_states> {{ {config_json['initial_state_x_n']} }};")
    else:
        lines.append("        x_n[ch].setZero();")
    if 'initial_state_v_n' in config_json:
        lines.append(f"        v_n[ch] = Eigen::Vector<T, num_nl_ports> {{ {config_json['initial_state_v_n']} }};")
    else:
        lines.append("        v_n[ch].setZero();")
    lines.append("    }")

    return lines

def generate_update_pots(netlist_info: NetlistInfo):
    u_fixed = []
    for el in netlist_info.elements:
        if el.element is Element.Voltage:
            if el.is_constant:
                u_fixed.append(el.name)
    u_fixed_str = ', '.join(u_fixed)

    lines = []
    lines.append(f"    Eigen::Matrix<T, num_pots, num_pots> Rv = Eigen::Matrix<T, num_pots, num_pots>::Zero();")
    for idx in range(netlist_info.num_pots):
        lines.append(f"    Rv ({idx}, {idx}) = pot_values[{idx}];")    
    lines.append("    Eigen::Matrix<T, num_pots, num_pots> Rv_Q_inv = (Rv + Q).inverse();")
    lines.append("")
    lines.append("    A_mat = A0_mat - (two_Z_Gx * (Ux * (Rv_Q_inv * Ux.transpose())));")
    lines.append("    Eigen::Matrix<T, num_states, num_voltages> B_mat = B0_mat - (two_Z_Gx * (Ux * (Rv_Q_inv * Uu.transpose())));")
    lines.append("    B_mat_var = B_mat.leftCols<num_voltages_variable>();")
    lines.append(f"    B_u_fix = B_mat.rightCols<num_voltages_constant>() * Eigen::Vector<T, num_voltages_constant> {{ {u_fixed_str} }};")
    lines.append("    C_mat = C0_mat - (two_Z_Gx * (Ux * (Rv_Q_inv * Un.transpose())));")
    lines.append("    D_mat = D0_mat - (Uo * (Rv_Q_inv * Ux.transpose()));")
    lines.append("    Eigen::Matrix<T, num_outputs, num_voltages> E_mat = E0_mat - (Uo * (Rv_Q_inv * Uu.transpose()));")
    lines.append("    E_mat_var = E_mat.leftCols<num_voltages_variable>();")
    lines.append(f"    E_u_fix = E_mat.rightCols<num_voltages_constant>() * Eigen::Vector<T, num_voltages_constant> {{ {u_fixed_str} }};")
    lines.append("    F_mat = F0_mat - (Uo * (Rv_Q_inv * Un.transpose()));")
    lines.append("    G_mat = G0_mat - (Un * (Rv_Q_inv * Ux.transpose()));")
    lines.append("    Eigen::Matrix<T, num_nl_ports, num_voltages> H_mat = H0_mat - (Un * (Rv_Q_inv * Uu.transpose()));")
    lines.append("    H_mat_var = H_mat.leftCols<num_voltages_variable>();")
    lines.append(f"    H_u_fix = H_mat.rightCols<num_voltages_constant>() * Eigen::Vector<T, num_voltages_constant> {{ {u_fixed_str} }};")
    lines.append("    K_mat = K0_mat - (Un * (Rv_Q_inv * Un.transpose()));")
    return lines

def generate_process_setup():
    lines = []
    lines.append("    Eigen::Vector<T, num_voltages_variable> u_n_var;")
    lines.append("    Eigen::Vector<T, num_nl_ports> p_n;")
    lines.append("    Eigen::Matrix<T, num_nl_ports, num_nl_ports> Jac = Eigen::Matrix<T, num_nl_ports, num_nl_ports>::Zero();")
    lines.append("    Eigen::Vector<T, num_nl_ports> i_n;")
    lines.append("    Eigen::Vector<T, num_nl_ports> F_min;")
    lines.append("    Eigen::Matrix<T, num_nl_ports, num_nl_ports> A_solve;")
    lines.append("    const Eigen::Matrix<T, num_nl_ports, num_nl_ports> eye = Eigen::Matrix<T, num_nl_ports, num_nl_ports>::Identity();")
    lines.append("    Eigen::Vector<T, num_nl_ports> delta_v;")
    lines.append("    Eigen::Vector<T, num_outputs> y_n;")
    lines.append("")

    return lines

def generate_process_sample(config_json, netlist_info: NetlistInfo, process_dtype):
    lines = []

    # TODO: what if there's more than one input, or it it's not index 0?
    lines.append("        u_n_var (0) = (T) sample;")
    lines.append("        p_n.noalias() = G_mat * x_n[ch] + H_mat_var * u_n_var + H_u_fix;")
    lines.append("")

    if netlist_info.num_nl_ports > 0:
        nl_elements = []
        for el in netlist_info.elements:
            if el.element in [Element.Diode, Element.DiodePair, Element.Transistor]:
                nl_elements.append(el)

        for el in nl_elements:
            if el.element is Element.DiodePair:
                lines.append("        const auto sinh_cosh = [] (T x)\n        {")
                lines.append("            // ref: https://en.wikipedia.org/wiki/Hyperbolic_functions#Definitions")
                lines.append("            // sinh = (e^(2x) - 1) / (2e^x), cosh = (e^(2x) + 1) / (2e^x)")
                lines.append("            // let B = e^x, then sinh = (B^2 - 1) / (2B), cosh = (B^2 + 1) / (2B)")
                lines.append("            // simplifying, we get: sinh = 0.5 (B - 1/B), cosh = 0.5 (B + 1/B)")
                lines.append("            auto B = std::exp (x);")
                lines.append("            auto Br = (T) 0.5 / B;")
                lines.append("            B *= (T) 0.5;")
                lines.append("            auto sinh = B - Br;")
                lines.append("            auto cosh = B + Br;")
                lines.append("            return std::make_pair (sinh, cosh);")
                lines.append("        };\n")
                break

        nl_idx = 0
        for el in nl_elements:
            if el.element is Element.Transistor:
                lines.append(f"        T exp_v{nl_idx+1}_v{nl_idx};")
                lines.append(f"        T exp_mv{nl_idx};")
                nl_idx += 2
            elif el.element is Element.Diode:
                lines.append(f"        T exp_v{nl_idx};")
                nl_idx += 1
            elif el.element is Element.DiodePair:
                lines.append(f"        T sinh_v{nl_idx};")
                lines.append(f"        T cosh_v{nl_idx};")
                nl_idx += 1

        lines.append("        const auto calc_currents = [&]\n        {")
        nl_idx = 0
        for el in nl_elements:
            if el.element is Element.Transistor:
                lines.append(f"            exp_v{nl_idx+1}_v{nl_idx} = std::exp ((v_n[ch]({nl_idx+1}) - v_n[ch]({nl_idx})) / Vt);")
                lines.append(f"            exp_mv{nl_idx} = std::exp (-v_n[ch]({nl_idx}) / Vt);")
                lines.append(f"            i_n ({nl_idx}) = Is_{el.name} * ((exp_v{nl_idx+1}_v{nl_idx} - (T) 1) / BetaF_{el.name} + (exp_mv{nl_idx} - (T) 1) / BetaR_{el.name});")
                lines.append(f"            i_n ({nl_idx+1}) = -Is_{el.name} * (-(exp_mv{nl_idx} - (T) 1) + AlphaF_{el.name} * (exp_v{nl_idx+1}_v{nl_idx} - (T) 1));")
                nl_idx += 2
            elif el.element is Element.Diode:
                lines.append(f"            exp_v{nl_idx} = std::exp (v_n[ch]({nl_idx}) / (n_{el.name} * Vt));")
                lines.append(f"            i_n ({nl_idx}) = -Is_{el.name} * (exp_v{nl_idx} - (T) 1);")
                nl_idx += 1
            elif el.element is Element.DiodePair:
                lines.append(f"            std::tie (sinh_v{nl_idx}, cosh_v{nl_idx}) = sinh_cosh (v_n[ch]({nl_idx}) / (n_{el.name} * Vt));")
                lines.append(f"            i_n ({nl_idx}) = (T) 2 * -Is_{el.name} * sinh_v{nl_idx};")
                nl_idx += 1
            lines.append("")
        lines.pop()
        lines.append("        };\n")

        lines.append("        const auto calc_jacobian = [&]\n        {")
        nl_idx = 0
        for el in nl_elements:
            if el.element is Element.Transistor:
                lines.append(f"            Jac ({nl_idx}, {nl_idx}) = (Is_{el.name} / Vt) * (-exp_v{nl_idx+1}_v{nl_idx} / BetaF_{el.name} - exp_mv{nl_idx} / BetaR_{el.name});")
                lines.append(f"            Jac ({nl_idx}, {nl_idx+1}) = (Is_{el.name} / Vt) * (exp_v{nl_idx+1}_v{nl_idx} / BetaF_{el.name});")
                lines.append(f"            Jac ({nl_idx+1}, {nl_idx}) = (Is_{el.name} / Vt) * (-exp_mv{nl_idx} + AlphaF_{el.name} * exp_v{nl_idx+1}_v{nl_idx});")
                lines.append(f"            Jac ({nl_idx+1}, {nl_idx+1}) = (Is_{el.name} / Vt) * (-AlphaF_{el.name} * exp_v{nl_idx+1}_v{nl_idx});")
                nl_idx += 2
            elif el.element is Element.Diode:
                lines.append(f"            Jac ({nl_idx}, {nl_idx}) = (-Is_{el.name} / (n_{el.name} * Vt)) * exp_v{nl_idx};")
                nl_idx += 1
            elif el.element is Element.DiodePair:
                lines.append(f"            Jac ({nl_idx}, {nl_idx}) = ((T) 2 * -Is_{el.name} / (n_{el.name} * Vt)) * cosh_v{nl_idx};")
                nl_idx += 1
            lines.append("")
        lines.pop()
        lines.append("        };\n")

        lines.append("        T delta;")
        lines.append("        int nIters = 0;")
        lines.append("        do\n        {")
        lines.append("            calc_currents();")
        lines.append("            calc_jacobian();\n")
        lines.append("            F_min.noalias() = p_n + K_mat * i_n - v_n[ch];")
        lines.append("            A_solve.noalias() = K_mat * Jac - eye;")
        lines.append("            delta_v.noalias() = A_solve.householderQr().solve (F_min);")
        lines.append("            v_n[ch] -= delta_v;")
        lines.append("            delta = delta_v.array().abs().sum();")
        lines.append(f"        }} while ({config_json['nr_exit_condition']});")
        lines.append("        calc_currents();")
        lines.append("")

    lines.append("        y_n.noalias() = D_mat * x_n[ch] + E_mat_var * u_n_var + E_u_fix + F_mat * i_n;")
    lines.append(f"        sample = ({process_dtype}) y_n (0);")
    lines.append("        x_n[ch] = A_mat * x_n[ch] + B_mat_var * u_n_var + B_u_fix + C_mat * i_n;")

    return lines


def generate_cpp(config_json, netlist_info: NetlistInfo, outputs):
    cpp_file = []

    cpp_file.append("/*")
    cpp_file.append(f" * This file was generated on {datetime.now()}")
    cpp_file.append(f" * using the command: `{' '.join(sys.argv)}`")
    cpp_file.append(" */")

    struct_name = config_json['struct_name']
    cpp_file.append(f"#include \"{struct_name}.h\"\n")

    cpp_file.append("namespace\n{")
    cpp_file.append("// START USER ENTRIES")
    for include_line in config_json['cpp_namespace_entries']:
        cpp_file.append(include_line)
    cpp_file.append("// END USER ENTRIES")
    cpp_file.append("} // namespace\n")

    if 'cpp_namespace' in config_json:
        cpp_file.append(f"namespace {config_json['cpp_namespace']} {{")

    # reset method
    cpp_file.append(f"void {struct_name}::reset (T fs)\n{{")
    cpp_file.append(generate_Gr(netlist_info))
    cpp_file.append(generate_Gx_Z(netlist_info))
    cpp_file.extend(generate_N_matrices(netlist_info, outputs))
    cpp_file.extend(generate_NDK_matrices(netlist_info))
    cpp_file.extend(generate_reset_states(config_json))
    cpp_file.append("}\n")

    # update method
    if netlist_info.num_pots > 0:
        cpp_file.append(f"void {struct_name}::update_pots (const std::array<T, num_pots>& pot_values)\n{{")
        cpp_file.extend(generate_update_pots(netlist_info))
        cpp_file.append("}\n")

    # process method
    process_dtype = config_json['process_data_type'] if 'process_data_type' in config_json else "T"
    cpp_file.append(f"void {struct_name}::process (std::span<{process_dtype}> channel_data, size_t ch) noexcept\n{{")
    cpp_file.extend(generate_process_setup())
    cpp_file.append("    for (auto& sample : channel_data)\n    {")
    cpp_file.extend(generate_process_sample(config_json, netlist_info, process_dtype))
    cpp_file.append("    }")
    cpp_file.append("}")

    if 'cpp_namespace' in config_json:
        cpp_file.append(f"}} // namespace {config_json['cpp_namespace']}")

    return cpp_file


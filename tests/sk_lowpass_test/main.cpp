#include <plt/matplotlibcpp.h>
#include "../freq_response_helpers.h"

namespace plt = matplotlibcpp;

#include "SKLowpassNDK.h"

namespace
{
constexpr double sample_rate = 48000.0;
constexpr double start_freq = 20.0;
constexpr double stop_freq = 20'000.0;
constexpr size_t N = 1 << 15;

std::array<double, 2> calc_R1_R2 (double fc, double Q)
{
    const auto C_base = std::sqrt (SKLowpassNDK::C1 * SKLowpassNDK::C2);
    const auto R_base = 1.0 / (2.0 * M_PI * C_base * fc);

    const auto m = (1.0 - std::sqrt (1.0 - 4.0 * Q * Q)) / (2.0 * Q);
    return { R_base * m, R_base / m };
}
} // namespace

int main()
{
    std::vector<double> input_data;
    input_data.resize (N);
    freq_helpers::generate_log_sweep (input_data, start_freq, stop_freq, sample_rate);

    std::vector<double> output_data (128, 0.0);

    SKLowpassNDK model;
    model.reset (sample_rate);
    model.update_pots (calc_R1_R2 (1000.0, 0.5));

    do
    {
        std::fill (output_data.begin(), output_data.end(), 0.0);
        model.process (output_data, 0);
    } while (std::abs (output_data[0]) > 1.0e-5);

    output_data = std::vector<double> { input_data.begin(), input_data.end() };
    model.process (output_data, 0);

    plt::semilogx (freq_helpers::fft_freqs (N / 2 + 1, sample_rate),
                   freq_helpers::compute_frequency_response (input_data, output_data));
    plt::grid (true);
    plt::xlim (start_freq, stop_freq);
    plt::save ("sk_lowpass.png");

    return 0;
}

#include <plt/matplotlibcpp.h>
#include "../freq_response_helpers.h"

namespace plt = matplotlibcpp;

#include "RCLowpassNDK.h"

namespace
{
constexpr double sample_rate = 48000.0;
constexpr double start_freq = 20.0;
constexpr double stop_freq = 20'000.0;
constexpr size_t N = 1 << 15;
} // namespace

int main()
{
    std::vector<double> input_data;
    input_data.resize (N);
    freq_helpers::generate_log_sweep (input_data, start_freq, stop_freq, sample_rate);

    std::vector<double> output_data (128, 0.0);

    RCLowpassNDK model;
    model.reset (sample_rate);
    model.update_pots ({ 1.0 / (2.0 * M_PI * RCLowpassNDK::C1 * 1000.0), 1.0 / (2.0 * M_PI * RCLowpassNDK::C2 * 1000.0) });

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
    plt::save ("rc_lowpass.png");

    return 0;
}

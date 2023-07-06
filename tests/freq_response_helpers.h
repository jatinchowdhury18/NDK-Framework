#pragma once

#include <cmath>
#include <complex>
#include <numeric>
#include <span>
#include <vector>

#include <kiss_fftr.h>

namespace freq_helpers
{
static void generate_log_sweep (std::span<double> sweep_buffer, double start_freq, double stop_freq, double sample_rate)
{
    const auto beta = (double) sweep_buffer.size() / std::log (stop_freq / start_freq);

    for (size_t i = 0; i < sweep_buffer.size(); i++)
    {
        double phase = 2.0 * (double) M_PI * beta * start_freq * (std::pow (stop_freq / start_freq, (double) i / (double) sweep_buffer.size()) - 1.0);

        sweep_buffer[i] = std::sin ((phase + (double) M_PI / 180.0) / sample_rate);
    }
}

double gain_to_db (double x)
{
    return 20.0 * std::log10 (std::max (x, 1.0e-5));
}

std::vector<double> compute_frequency_response (std::span<const double> sweep_buffer, std::span<const double> filter_buffer)
{
    const auto fft_data_size = sweep_buffer.size();
    auto fft = kiss_fftr_alloc ((int) fft_data_size, 0, nullptr, nullptr);

    std::vector<std::complex<double>> sweep_fft;
    sweep_fft.resize (sweep_buffer.size() / 2 + 1);
    kiss_fftr (fft, sweep_buffer.data(), reinterpret_cast<kiss_fft_cpx*> (sweep_fft.data()));

    std::vector<std::complex<double>> filter_fft;
    filter_fft.resize (filter_buffer.size() / 2 + 1);
    kiss_fftr (fft, filter_buffer.data(), reinterpret_cast<kiss_fft_cpx*> (filter_fft.data()));

    kiss_fftr_free (fft);

    std::vector<double> mag_response_dB (sweep_buffer.size() / 2 + 1, 0.0);
    for (size_t i = 0; i < mag_response_dB.size(); ++i)
        mag_response_dB[i] = gain_to_db (std::abs (filter_fft[i] / sweep_fft[i]));

    return mag_response_dB;
}

static std::vector<double> fft_freqs (size_t num, double sample_rate)
{
    auto val = 0.5 * sample_rate / (double) num;

    std::vector<double> results (num, 0.0);
    std::iota (results.begin(), results.end(), 0.0);
    std::transform (results.begin(), results.end(), results.begin(), [val] (auto x)
                    { return x * val; });

    return results;
}
} // namespace freq_helpers

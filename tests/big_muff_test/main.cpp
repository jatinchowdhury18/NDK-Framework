#include <numeric>

#include <plt/matplotlibcpp.h>
namespace plt = matplotlibcpp;

#include "BigMuff2D.h"
#include "BigMuffDP.h"

namespace
{
    constexpr double fs = 48000.0;
    constexpr double fc = 10.0;
    constexpr size_t N = 240000;
}

int main()
{
    std::vector<double> input_data (N, 0.0);
    for (size_t n = 0; n < N; ++n)
        input_data[n] = 1.0 * std::sin (2.0 * M_PI * (float) n * fc / fs) + 1.0;

    std::vector<double> output_data_2d { input_data.begin(), input_data.end() };
    BigMuff2D model_2d;
    model_2d.reset (fs);
    model_2d.process (output_data_2d, 0);

    std::vector<double> output_data_dp { input_data.begin(), input_data.end() };
    BigMuffDP model_dp;
    model_dp.reset (fs);
    model_dp.process (output_data_dp, 0);

    plt::plot (input_data);
    plt::plot (output_data_2d);
    plt::plot (output_data_dp, "--");
    plt::save ("big_muff.png");
}

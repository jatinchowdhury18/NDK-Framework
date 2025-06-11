#include <numeric>

#include <plt/matplotlibcpp.h>
namespace plt = matplotlibcpp;

#include "CommonCathodeNDK.h"

namespace
{
    constexpr double fs = 96000.0;
    constexpr double fc = 1000.0;
    constexpr size_t N = 4800;
}

int main()
{
    std::vector<double> input_data (N, 0.0);
    for (size_t n = 0; n < N; ++n)
        input_data[n] = 0.5 * std::sin (2.0 * M_PI * (float) n * fc / fs) + 0.0;

    std::vector<double> output_data { input_data.begin(), input_data.end() };
    CommonCathodeNDK model;
    model.reset (fs);
    model.process (output_data, 0);

    plt::plot (input_data);
    plt::plot (output_data);
//    plt::save ("common_cathode.png");
    plt::show();
}

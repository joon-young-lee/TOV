#include <iostream>
#include "CubicSpline.hpp"

int main() {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 1.0, 4.0, 9.0, 16.0};

    CubicSpline spline(x, y);

    std::vector<double> test_points = {0.5, 1.5, 2.5, 3.5};

    for (double point : test_points) {
        std::cout << "Spline(" << point << ") = " << spline(point) << std::endl;
    }

    return 0;
}

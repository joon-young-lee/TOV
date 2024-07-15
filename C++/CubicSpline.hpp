#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include <vector>

class CubicSpline {
public:
    CubicSpline(const std::vector<double>& x, const std::vector<double>& y);
    double operator()(double x) const;

private:
    void computeCoefficients();

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
};

#endif // CUBIC_SPLINE_HPP

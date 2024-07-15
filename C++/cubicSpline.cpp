#include "CubicSpline.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>

CubicSpline::CubicSpline(const std::vector<double> &x, const std::vector<double> &y)
    : x_(x), y_(y), a_(x.size()), b_(x.size()), c_(x.size()), d_(x.size())
{
    if (x.size() != y.size() || x.size() < 2)
    {
        throw std::invalid_argument("Input vectors must have the same size and contain at least two points.");
    }

    computeCoefficients();
}

void CubicSpline::computeCoefficients()
{
    size_t n = x_.size() - 1;
    std::vector<double> h(n), alpha(n), l(n + 1), mu(n), z(n + 1);

    for (size_t i = 0; i < n; ++i)
    {
        h[i] = x_[i + 1] - x_[i];
    }

    for (size_t i = 1; i < n; ++i)
    {
        alpha[i] = (3 / h[i]) * (y_[i + 1] - y_[i]) - (3 / h[i - 1]) * (y_[i] - y_[i - 1]);
    }

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for (size_t i = 1; i < n; ++i)
    {
        l[i] = 2 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1.0;
    z[n] = 0.0;
    c_[n] = 0.0;

    for (size_t j = n - 1; j != static_cast<size_t>(-1); --j)
    {
        c_[j] = z[j] - mu[j] * c_[j + 1];
        b_[j] = (y_[j + 1] - y_[j]) / h[j] - h[j] * (c_[j + 1] + 2 * c_[j]) / 3;
        d_[j] = (c_[j + 1] - c_[j]) / (3 * h[j]);
        a_[j] = y_[j];
    }
}

double CubicSpline::operator()(double x) const
{
    if (x < x_.front() || x > x_.back())
    {
        throw std::out_of_range("Extrapolation not supported");
    }

    size_t n = x_.size() - 1;
    size_t i = 0;
    for (size_t j = 0; j < n; ++j)
    {
        if (x_[j] <= x && x < x_[j + 1])
        {
            i = j;
            break;
        }
    }

    double dx = x - x_[i];
    return a_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
}
#include <iostream>
#include <vector>
#include "CubicSpline.hpp"

// Constants (define them appropriately as per your requirement)
const double G = 6.67430e-8; // Gravitational constant in cgs
const double c = 2.998e10;   // Speed of light in cm/s
const double pi = 3.141592653589793;

double MeV_fm_to_cgs = 1.60218e33; // Example conversion factor, adjust as needed

// TOV solver function
void TOV(const CubicSpline& EoS, double p0, double del_r, int maxit) {
    auto dp_dr = [&EoS](double r, double p, double m) -> double {
        double e = EoS(p);
        return -G * e * m * (1 + p / e) * (1 + 4 * pi * p * std::pow(r, 3) / (m * c * c)) /
               (std::pow(c * r, 2) * (1 - 2 * G * m / (c * c * r)));
    };

    auto dm_dr = [&EoS](double r, double p, double m) -> double {
        double e = EoS(p);
        return 4 * pi * e * std::pow(r, 2) / (c * c);
    };

    double r = 1.e-10;
    double p = p0;
    double m = 1.e-30;
    int i = 0;

    while (i <= maxit && p >= 5.38125462909e15) {
        double k1_p = dp_dr(r, p, m);
        double k1_m = dm_dr(r, p, m);

        double r2 = r + del_r / 2;
        double p2 = p + del_r * k1_p / 2;
        double m2 = m + del_r * k1_m / 2;
        double k2_p = dp_dr(r2, p2, m2);
        double k2_m = dm_dr(r2, p2, m2);

        double r3 = r2;
        double p3 = p + del_r * k2_p / 2;
        double m3 = m + del_r * k2_m / 2;
        double k3_p = dp_dr(r3, p3, m3);
        double k3_m = dm_dr(r3, p3, m3);

        double r4 = r + del_r;
        double p4 = p + del_r * k3_p;
        double m4 = m + del_r * k3_m;
        double k4_p = dp_dr(r4, p4, m4);
        double k4_m = dm_dr(r4, p4, m4);

        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r;
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r;
        r += del_r;
        i += 1;

        std::cout << "Mass: " << m << ", Pressure: " << p << std::endl;
        std::cout << "_______________________________" << std::endl;
    }
}

// Main function
int main() {
    // Load your data (example data)
    std::vector<double> p_arr = {1e34, 2e34, 3e34, 4e34}; // Example pressure data in cgs
    std::vector<double> e_arr = {1e15, 2e15, 3e15, 4e15}; // Example energy density data in cgs

    // Convert units if necessary (assuming the data needs to be converted from MeV/fm^3 to cgs)
    for (auto& p : p_arr) p *= MeV_fm_to_cgs;
    for (auto& e : e_arr) e *= MeV_fm_to_cgs;

    // Create the cubic spline interpolation for the equation of state
    CubicSpline EoS(p_arr, e_arr);

    // Set initial pressure, step size, and max iterations
    double p0 = p_arr.front();
    double del_r = 1.0; // Example step size in cm
    int maxit = 1000; // Example maximum iterations

    // Solve the TOV equations
    TOV(EoS, p0, del_r, maxit);

    return 0;
}

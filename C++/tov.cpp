#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// Constants in CGS units
const double G = 6.67430e-8;  // Gravitational constant (cm^3/g/s^2)
const double c = 2.99792458e10;  // Speed of light (cm/s)
const double c2 = c * c;  // Speed of light squared (cm^2/s^2)
const double M0 = 1.9884e33;  // Solar mass (g)
const double MeV_fm_to_cgs = 1.602179e-6;  // Conversion factor from MeV/fm^3 to CGS units (erg/cm^3)
const double erg_cm_to_MeV_fm = 6.24151e5;  // Conversion factor from erg/cm^3 to MeV/fm^3

// Cubic spline interpolation
class CubicSpline {
public:
    CubicSpline(const std::vector<double>& x, const std::vector<double>& y, bool extrapolate = false) {
        // Implement cubic spline interpolation
    }

    double operator()(double x) const {
        // Implement cubic spline evaluation
    }
};

// Equation of state
std::pair<CubicSpline, std::vector<double>> eos(const std::string& file_name) {
    std::ifstream file(file_name);
    std::vector<double> e_arr, p_arr;
    double e, p;
    while (file >> e >> p) {
        e_arr.push_back(e * MeV_fm_to_cgs);
        p_arr.push_back(p * MeV_fm_to_cgs);
    }
    std::reverse(e_arr.begin(), e_arr.end());
    std::reverse(p_arr.begin(), p_arr.end());
    return {CubicSpline(p_arr, e_arr, false), p_arr};
}

// Tolman-Oppenheimer-Volkoff (TOV) equation
std::pair<double, double> tov(const CubicSpline& eos, double p0, double del_r, int maxit) {
    double r = 1e-10, p = p0, m = 1e-30;
    int i = 0, cut = 6000;

    auto dp_dr = [&](double r, double p, double m) {
        return -G * eos(p) * m * (1 + p / eos(p)) * (1 + 4 * M_PI * p * r * r * r / (m * c2)) / (1 - 2 * G * m / (c2 * r)) / (c * r * r);
    };

    auto dm_dr = [&](double r, double p, double m) {
        return 4 * M_PI * eos(p) * r * r / c2;
    };

    for (; i < cut; ++i) {
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

        double r4 = r2 + del_r / 2;
        double p4 = p + del_r * k3_p;
        double m4 = m + del_r * k3_m;
        double k4_p = dp_dr(r4, p4, m4);
        double k4_m = dm_dr(r4, p4, m4);

        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r;
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r;
        r += del_r;
    }

    while (i <= maxit && p >= 5.38125462909e23) {
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

        double r4 = r3 + del_r / 2;
        double p4 = p + del_r * k3_p;
        double m4 = m + del_r * k3_m;
        double k4_p = dp_dr(r4, p4, m4);
        double k4_m = dm_dr(r4, p4, m4);

        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r;
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r;
        r += del_r;
        ++i;
    }

    return {m / M0, r / 1e5};
}

int main() {
    double del_r = 100.0;
    int maxit = 1e5;

    for (int j = 50; j <= 100; ++j) {
        std::string file_name = "Archive/eft_pnm32_" + std::to_string(j).substr(0, 6) + "_ldm_eos.dat";
        auto [eos, p_arr] = eos(file_name);
        int step = 1000;
        std::vector<double> M(step), R(step), P(step);

        double p_i = 1.602179e+33;  // 1.0 MeV/fm^3
        double p_f = 1.602179e+36;  // 1000 MeV/fm^3
        for (int l = 1; l <= step; ++l) {
            P[l - 1] = p_i + (l - 1) * (p_f - p_i) / step;
            std::tie(M[l - 1], R[l - 1]) = tov(eos, P[l - 1], del_r, maxit);
            P[l - 1] *= erg_cm_to_MeV_fm;  // MeV/fm^3
            std::cout << "-------------------------------------------" << std::endl;
            std::cout << "Archive/eft_pnm32_" << std::setw(6) << std::setfill('0') << j << "_ldm_eos.dat: " << l / 10.0 << "%" << std::endl;
            std::cout << "Initial pressure: " << P[l - 1] << " MeV/fm^3" << std::endl;
            std::cout << M[l - 1] << " M_0" << std::endl;
            std::cout << R[l - 1] << " km" << std::endl;
            std::cout << "-------------------------------------------" << std::endl;
        }

        std::cout << "Initial Pressure " << p_i << " ~ " << p_f << std::endl;
        std::ofstream out("Result_2/eft_pnm32_" + std::to_string(j).substr(0, 6) + "_ldm_result.dat");
        for (int i = 0; i < step; ++i) {
            out << std::scientific << std::setprecision(6) << P[i] << " " << M[i] << " " << R[i] << std::endl;
        }
    }

    return 0;
}

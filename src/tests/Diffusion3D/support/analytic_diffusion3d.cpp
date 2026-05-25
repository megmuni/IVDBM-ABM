#include "analytic_diffusion3d.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

} // namespace

std::vector<double> analytic_gaussian_impulse_3d(int nx, int ny, int nz,
                                                   double h, double D,
                                                   double t,
                                                   int cx, int cy, int cz,
                                                   double mass)
{
    if (t <= 0.0)
        throw std::invalid_argument("analytic_gaussian_impulse_3d: t must be > 0");
    if (D <= 0.0)
        throw std::invalid_argument("analytic_gaussian_impulse_3d: D must be > 0");

    const int n = nx * ny * nz;
    std::vector<double> out(static_cast<size_t>(n), 0.0);

    const double denom = std::pow(4.0 * M_PI * D * t, 1.5);
    const double coeff = mass / denom;
    const double inv4Dt = 1.0 / (4.0 * D * t);
    const double h2 = h * h;

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const int ix = x - cx;
                const int iy = y - cy;
                const int iz = z - cz;
                const double r2 = h2 * static_cast<double>(ix * ix + iy * iy + iz * iz);
                out[static_cast<size_t>(x + nx * (y + ny * z))] =
                    coeff * std::exp(-r2 * inv4Dt);
            }
        }
    }

    return out;
}

double rel_l2_error(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("rel_l2_error: size mismatch");

    long double num = 0.0L;
    long double den = 0.0L;
    for (size_t i = 0; i < a.size(); ++i)
    {
        const long double d = static_cast<long double>(a[i]) - static_cast<long double>(b[i]);
        num += d * d;
        den += static_cast<long double>(b[i]) * static_cast<long double>(b[i]);
    }

    if (den == 0.0L)
        return (num == 0.0L) ? 0.0 : std::numeric_limits<double>::infinity();

    return std::sqrt(static_cast<double>(num / den));
}

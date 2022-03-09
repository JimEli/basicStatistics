#ifndef NORMAL_H
#define NORMAL_H

#include "common.h"
#include "erf.h"

/*
  The Normal Distribution
  Distribution and quantile functions for the normal distribution.
    dnorm returns the probability density function.
    pnorm gives the distribution function.
    qnorm gives the quantile function.
    (pnorm and qnorm are inverses of each other)
  If mu or sdev are not specified they assume the default values of 0 and 1, respectively.
  The normal distribution has density:
     f(x) = 1 / (sqrt(2 * PI) sigma) e ^ -((x - mu) ^ 2 / (2 sigma ^ 2))
  Usage:
     dnorm(x, mean = 0, sd = 1)
     pnorm(q, mean = 0, sd = 1)
     qnorm(p, mean = 0, sd = 1)
  Arguments:
     x, q = quantile
     p = probability
     n = number of observations, if length(n)>1, the length is taken to be the number required
     mu = mean
     sdev = standard deviation
*/

// Returns probability/percent/proportion/area under curve of normal distribution [0 to x].
// If n > 30 can use normal distribution.
// If n <= 30, then original population must be normal distribution.
double pNorm(const double n, const double mu = 0., const double sdev = 1.)
{
    double z = ((n - mu) / sdev);
    return _erfc(-z / sqrt(2.)) / 2.0;
}

// Returns x of given probability/percent/proportion/area under curve of normal distribution.
// If n*p > 5 and n *(1 - p) > 5, then can use normal distribution.
double qNorm(const double p, const double mu = 0., const double sdev = 1.)
{
    double z = _erf(p);
    return mu + (z * sdev);
}

// Compute the density of the normal distribution.
double dNorm(const double x, const double mu = 0., const double sigma = 1.)
{
  //return (1. / (sigma * sqrt(2. * 3.14159))) * pow(2.71828, (pow(x - mu, 2.) / (2. * sigma * sigma)));

#ifdef IEEE_754
  if (isnan(x) || isnan(mu) || isnan(sigma))
    return x + mu + sigma;
#endif
  if (sigma < 0)
    return NAN;
    
  if (!isfinite(sigma))
    return 0.;
    
  if (!isfinite(x) && mu == x)
    return NAN;

  if (sigma == 0)
    return (x == mu) ? ML_POSINF : 0.;
    
  double x_ = (x - mu) / sigma;

  if (!isfinite(x_))
    return 0.;

  x_ = fabs(x_);
  if (x_ >= 2 * sqrt(DBL_MAX))
    return 0.;

  return M_1_SQRT_2PI * exp(-0.5 * x_ * x_) / sigma;
}

// Cumulative normal distribution.
double pNormCDF(double x)
{
    // return 0.5 * erfc(-value * M_SQRT1_2);
    constexpr double a1 = 0.254829592;
    constexpr double a2 = -0.284496736;
    constexpr double a3 = 1.421413741;
    constexpr double a4 = -1.453152027;
    constexpr double a5 = 1.061405429;
    constexpr double p = 0.3275911;

    // Save sign of x.
    int sign = 1;
    if (x < 0)
        sign = -1;

    x = fabs(x) / sqrt(2.);

    double t = 1. / (1. + p * x);
    double y = 1. - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1. + sign * y);
}

// Compute the quantile function for the normal distribution.
double qNormCDF(double p, double mu, double sigma)
{
    if (p < 0. || p > 1.)
    {
        std::cout << "The probality p must be bigger than 0 and smaller than 1\n";
        exit(1);
    }

    if (sigma < 0.)
    {
        std::cout << "The standard deviation sigma must be positive\n";
        exit(1);
    }

    if (p == 0.)
        return -ML_NEGINF;

    if (p == 1.)
        return ML_POSINF;

    if (sigma == 0.)
        return mu;

    double q, r, val;

    q = p - 0.5;

    if (fabs(q) <= .425)
    {
        r = .180625 - q * q;
        val = q * (((((((r * 2509.0809287301226727 +
            33430.575583588128105) * r + 67265.770927008700853) * r +
            45921.953931549871457) * r + 13731.693765509461125) * r +
            1971.5909503065514427) * r + 133.14166789178437745) * r +
            3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1);
    }
    else
    {
        if (q > 0)
            r = 1 - p;
        else
            r = p;

        r = sqrt(-log(r));

        if (r <= 5)
        {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                .0227238449892691845833) * r + .24178072517745061177) *
                r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
                r + 4.6303378461565452959) * r +
                1.42343711074968357734)
                / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1);
        }
        else
        {
            r += -5;
            val = (((((((r * 2.01033439929228813265e-7 +
                2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
                r + .29656057182850489123) * r +
                1.7848265399172913358) * r + 5.4637849111641143699) *
                r + 6.6579046435011037772)
                / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1);
        }

        if (q < 0.0)
        {
            val = -val;
        }
    }

    return mu + sigma * val;
}

#endif

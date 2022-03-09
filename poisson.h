#ifndef POISSON_H
#define POISSON_H

#include "common.h"

//double dPois(const double k, double lambda) { return (pow(lambda, k) * exp(-lambda)) / std::tgamma(k + 1); }
//double pPois(const double k, const double lambda) { return exp(k * log(lambda) - lgamma(k + 1.0) - lambda); }

//
// dPois
//

// Poisson distribution probability mass function (PMF).
//   Poisson probability function with mean lambda can be calculated for any value of x.
// Arguments:
//   x=axis values(x = 0, 1, 2, ...)
//   lambda=mean number of events that occur on the interval
double dPois(double x, double lambda)
{
#ifdef IEEE_754
    if (isnan(x) || isnan(lambda))
        return x + lambda;
#endif

    if (lambda < 0.)
        return NAN;

    if (x < 0 || !isfinite(x))
        return 0.;

    x = round(x);

    return dpois_raw(x, lambda);
}

//
// pPois
//

// Compute the log of a sum from logs of terms, i.e.,
double logspace_add(double logx, double logy) { return fmax2(logx, logy) + log1p(exp(-fabs(logx - logy))); }
// Compute the log of a difference from logs of terms, i.e.,
double logspace_sub(double logx, double logy) { return logx + ((logy - logx) > -M_LN2 ? log(-expm1(logy - logx)) : log1p(-exp(logy - logx))); }
// Compute the log of a sum from logs of terms, i.e.,
double logspace_sum(const double* logx, int n)
{
    if (n == 0)
        return ML_NEGINF; // = log( sum(<empty>) )

    if (n == 1)
        return logx[0];

    if (n == 2)
        return logspace_add(logx[0], logx[1]);

    int i;
    double Mx = logx[0];

    for (i = 1; i < n; i++)
        if (Mx < logx[i]) Mx = logx[i];

    long double s = (long double)0.;

    for (i = 0; i < n; i++)
        s += expl(logx[i] - Mx);

    return Mx + (double)logl(s);
}

// Compute the following ratio with higher accuracy that would be had
// from doing it directly.
//     dnorm (x, 0, 1, FALSE)
//     ----------------------------------
//     pnorm (x, 0, 1, lower_tail, FALSE)
static double dpnorm(double x, int lower_tail, double lp)
{
    if (x < 0)
    {
        x = -x;
        lower_tail = !lower_tail;
    }

    if (x > 10 && !lower_tail)
    {
        double term = 1 / x;
        double sum = term;
        double x2 = x * x;
        double i = 1;

        do {
            term *= -i / x2;
            sum += term;
            i += 2;
        } while (fabs(term) > DBL_EPSILON * sum);

        return 1 / sum;
    }
    else
    {
        double d;

        if (!lower_tail)
            d = 1. - dNorm(x, 0., 1.);
        else
            d = 1. - dNorm(x, 0., 1.);

        return d / exp(lp);
    }
}

// Poisson cumulative distribution function (CDF).
//   The probability of a variable X following a Poisson distribution 
//   taking values equal or lower than x can be calculated.
// Arguments:
//   q=quantile
//   lambda=mean
double pPois(double x, double lambda)
{
#ifdef IEEE_754
    if (isnan(x) || isnan(lambda))
        return x + lambda;
#endif

    if (lambda < 0.)
        return NAN;

    if (x < 0)
        return 0.;

    if (lambda == 0.)
        return 1.;

    if (!isfinite(x))
        return 1.;

    x = floor(x + 1e-7);

    int lower_tail = false;
    return pgamma(lambda, x + 1, 1., lower_tail);
}

// The quantile function of the Poisson distribution.
static double doPoisSearch(double y, double* z, double p, double lambda, double incr)
{
    if (*z >= p)
    {
        // search to left 
        for (;;)
        {
            if (y == 0 || (*z = pPois(y - incr, lambda)) < p)
                return y;
            y = fmax2(0, y - incr);
        }
    }
    else
    {
        // search to right 
        for (;;)
        {
            y = y + incr;
            if ((*z = pPois(y, lambda)) >= p)
                return y;
        }
    }
}

double qPois(double p, double lambda)
{
    double mu, sigma, gamma, z, y;

#ifdef IEEE_754
    if (isnan(p) || isnan(lambda))
        return p + lambda;
#endif

    if (!isfinite(lambda))
        return NAN;

    if (lambda < 0)
        return NAN;

    if (lambda == 0)
        return 0;

    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return 0;
    if (p == 1)
        return ML_POSINF;

    mu = lambda;
    sigma = sqrt(lambda);
    gamma = 1.0 / sigma;

    if (p + 1.01 * DBL_EPSILON >= 1.)
        return ML_POSINF;

    z = qNorm(p, 0., 1.);
    y = round(mu + sigma * (z + gamma * (z * z - 1) / 6));
    z = pPois(y, lambda);

    p *= 1 - 64 * DBL_EPSILON;

    if (lambda < 1e5)
        return doPoisSearch(y, &z, p, lambda, 1);

    double incr = floor(y * 0.001), oldincr;

    do
    {
        oldincr = incr;
        y = doPoisSearch(y, &z, p, lambda, incr);
        incr = fmax2(1, floor(incr / 100));
    } while (oldincr > 1 && incr > lambda * 1e-15);

    return y;
}

#endif

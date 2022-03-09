#ifndef BINOMIAL_H
#define BINOMIAL_H

#include "common.h"

/*
   The Binomial Distribution

   Density, quantile and distribution function for the binomial distribution, 
     dbinom gives the density
     pbinom gives the distribution function
     qbinom gives the quantile function
   This is conventionally interpreted as the number of ‘successes’ in size trials. 
   
   The binomial distribution with size = n and prob = p has density:
     p(x) = choose(n, x) p^x (1-p)^(n-x)
   for x = 0, …, n. Note that binomial coefficients can be computed by choose in R.
   The quantile is defined as the smallest value x such that F(x) ≥ p, where F is the distribution function.
   
   Usage:
     dbinom(x, size, prob, log = FALSE)
     pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
   
   Arguments:
     x, q=quantiles
     p=probability
     n=number of observations. If length(n) > 1, the length is taken to be the number required
     k=number of trials (zero or more).
*/

// Binomial PMF(p, n,k) = n!/(k!*(n-k)!) p^k (1-p)^(n-k)
static double logChoose(const unsigned n, const unsigned k)
{
    return std::lgamma(double(n + 1)) - std::lgamma(double(k + 1)) - std::lgamma(double(n - k + 1));
}

// Binomial PMF(p, n,k) = n!/(k!*(n-k)!) p^k (1-p)^(n-k)
double dBinom(const unsigned k, const unsigned n, const double p)
{
    //double nCk = (factorial((unsigned)n) / (factorial((unsigned)k) * factorial((unsigned)(n - k))));
    //return (nCk * pow(p, k) * pow(1.0 - p, n - k));
    double lgr = logChoose(n, k) + double(k) * std::log(p) + double(n - k) * std::log(1 - p);

    return std::exp(lgr);
}

// Binomial CDF.
double pBinom(const unsigned k, const unsigned n, const double p)
{
    double cdf = 0.;
    double b = 0.;

    for (unsigned _k = 1; _k <= (unsigned)k; _k++)
    {
        double log_pmf_k = 0.;

        b += +log(n - _k + 1.) - log(_k);
        log_pmf_k = b + _k * log(p) + (n - _k) * log(1. - p);
        cdf += exp(log_pmf_k);
    }

    return cdf;
}

static double doBinomSearch(double y, double* z, double p, double n, double pr, double incr)
{
    if (*z >= p) 
    {
        // Search to left.
        for (;;) 
        {
            double newz;

            if (y == 0 || (newz = pBinom((unsigned)(y - incr), (unsigned)n, pr)) < p)
                return y;
            y = fmax2(0, y - incr);
            *z = newz;
        }
    } 
    else 
    {		
        // Search to right.
        for (;;) 
        {
            y = fmin2(y + incr, n);
            if (y == n || (*z = pBinom((unsigned)y, (unsigned)n, pr)) >= p)
                return y;
        }
    }
}

// The quantile function of the binomial distribution.
double qBinom(double p, double n, double pr)
{
    double q, mu, sigma, gamma, z, y;

#ifdef IEEE_754
    if (isnan(p) || isnan(n) || isnan(pr))
        return p + n + pr;
#endif

    if (!isfinite(n) || !isfinite(pr))
        return NAN;

    if (!isfinite(p))
        return NAN;

    if (n != floor(n + 0.5))
        return NAN;

    if (pr < 0 || pr > 1 || n < 0)
        return NAN;

    if (p < 0 || p > 1)
      return NAN;
    if (p == 0)
      return 0;
    if (p == 1)
      return n;
  
    if (pr == 0. || n == 0)
        return 0.;

    q = 1 - pr;
    if (q == 0.)
        return n;

    mu = n * pr;
    sigma = sqrt(n * pr * q);
    gamma = (q - pr) / sigma;

    if (p + 1.01 * DBL_EPSILON >= 1.)
        return n;

    z = qNorm(p, 0., 1.);
    y = floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);

    if (y > n) 
        y = n;

    z = pBinom((unsigned)y, (unsigned)n, pr);

    p *= 1 - 64 * DBL_EPSILON;

    if (n < 1e5)
        return doBinomSearch(y, &z, p, n, pr, 1);

    double incr = floor(n * 0.001), oldincr;

    do 
    {
      oldincr = incr;
      y = doBinomSearch(y, &z, p, n, pr, incr);
      incr = fmax2(1, floor(incr / 100));
    } while (oldincr > 1 && incr > n * 1e-15);

    return y;
}

#endif

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

/*
static unsigned factorial(const unsigned n)
{
    //return std::tgamma(n + 1);
    unsigned ret = 1;
    for (unsigned i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}
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
        // search to the left 
        for (;;) 
        {
            double newz;

            //if (y == 0 || (newz = pbinom(y - incr, n, pr, true, false)) < p)
            if (y == 0 || (newz = pBinom((unsigned)(y - incr), (unsigned)n, pr)) < p)
                return y;
            y = fmax2(0, y - incr);
            *z = newz;
        }
    } 
    else 
    {		
        // search to the right 
        for (;;) 
        {
            y = fmin2(y + incr, n);
            //if (y == n || (*z = pbinom(y, n, pr, false)) >= p)
            if (y == n || (*z = pBinom((unsigned)y, (unsigned)n, pr)) >= p)
                return y;
        }
    }
}

// The quantile function of the binomial distribution.
// Uses the Cornish-Fisher Expansion to include a skewness correction to a normal approximation. This gives an
// initial value which never seems to be off by more than 1 or 2. A search is then conducted of values close to
// this initial start point.
double qBinom(double p, double n, double pr)
{
    double q, mu, sigma, gamma, z, y;

#ifdef IEEE_754
    if (isnan(p) || isnan(n) || isnan(pr))
        return p + n + pr;
#endif

    if (!isfinite(n) || !isfinite(pr))
        return NAN;

    // if log_p is true, p = -Inf is a legitimate value 
    //if (!isfinite(p) && !log_p) return NAN;
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
        return n; // covers the full range of the distribution 

    mu = n * pr;
    sigma = sqrt(n * pr * q);
    gamma = (q - pr) / sigma;

    // temporary hack --- FIXME --- 
    if (p + 1.01 * DBL_EPSILON >= 1.)
        return n;

    // y := approx.value (Cornish-Fisher expansion): 
    //z = qnorm(p, 0., 1., true, false);
    z = qNorm(p, 0., 1.);
    y = floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);

    if (y > n) // way off 
        y = n;

    //z = pbinom(y, n, pr, true, false);
    z = pBinom((unsigned)y, (unsigned)n, pr);

    // fuzz to ensure left continuity: 
    p *= 1 - 64 * DBL_EPSILON;

    if (n < 1e5)
        return doBinomSearch(y, &z, p, n, pr, 1);

    // Otherwise be a bit cleverer in the search 
    double incr = floor(n * 0.001), oldincr;

    do 
    {
      oldincr = incr;
      y = doBinomSearch(y, &z, p, n, pr, incr);
      incr = fmax2(1, floor(incr / 100));
    } while (oldincr > 1 && incr > n * 1e-15);

    return y;
}

#if 0

// Saddle point algorithm.

// NTYPE is the type used for the n and x arguments.
// For 32-bit integers, the maximum n is 2^31-1=2147483647.
// If larger n is required, NTYPE must be double.
typedef int NTYPE;

#define PI2 6.283185307179586476925286
#define S0 0.083333333333333333333        // 1/12 
#define S1 0.00277777777777777777778      // 1/360 
#define S2 0.00079365079365079365079365   // 1/1260 
#define S3 0.000595238095238095238095238  // 1/1680 
#define S4 0.0008417508417508417508417508 // 1/1188 

static double sfe[16] = 
{
  0, 0.081061466795327258219670264,
  0.041340695955409294093822081, 0.0276779256849983391487892927,
  0.020790672103765093111522771, 0.0166446911898211921631948653,
  0.013876128823070747998745727, 0.0118967099458917700950557241,
  0.010411265261972096497478567, 0.0092554621827127329177286366,
  0.008330563433362871256469318, 0.0075736754879518407949720242,
  0.006942840107209529865664152, 0.0064089941880042070684396310,
  0.005951370112758847735624416, 0.0055547335519628013710386899
};

// stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) 
double stirlerr(NTYPE n)
{ 
    double nn;

    if (n < 16) 
        return (sfe[(int)n]);
    nn = (double)n;
    nn = nn * nn;
    if (n > 500) 
        return ((S0 - S1 / nn) / n);
    if (n > 80) 
        return ((S0 - (S1 - S2 / nn) / nn) / n);
    if (n > 35) 
        return ((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
    return ((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
}

// Evaluate the deviance term bd0(x,np) = x log(x/np) + np - x
double bd0(NTYPE x, double np)
{ 
    double ej, s, s1, v;
    int j;

    if (fabs(x - np) < 0.1 * (x + np))
    {
        s = (x - np) * (x - np) / (x + np);
        v = (x - np) / (x + np);
        ej = 2 * x * v;
        for (j = 1; ; j++)
        {
            ej *= v * v;
            s1 = s + ej / (2 * j + 1);
            if (s1 == s) 
                return (s1);
            s = s1;
        }
    }

    return (x * log(x / np) + np - x);
}

double dbinom(NTYPE x, NTYPE n, double p)
{ 
    double lc;

    if (p == 0.0) 
        return ((x == 0) ? 1.0 : 0.0);
    if (p == 1.0) 
        return ((x == n) ? 1.0 : 0.0);
    if (x == 0) 
        return (exp(n * log(1 - p)));
    if (x == n) 
        return (exp(n * log(p)));
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * (1.0 - p));
    return(exp(lc) * sqrt(n / (PI2 * x * (n - x))));
}

double dpois(NTYPE x, double lb)
{ 
    if (lb == 0) 
        return ((x == 0) ? 1.0 : 0.0);
    if (x == 0) 
        return (exp(-lb));
    return (exp(-stirlerr(x) - bd0(x, lb)) / sqrt(PI2 * x));
}

// Mutiplication algorithm.
double dbinom_mult(int x, int n, double p;)
{ 
    double f;
    int j0, j1, j2;

    if (2 * x > n) 
        return (dbinom_mult(n - x, n, 1 - p));
    j0 = j1 = j2 = 0;
    f = 1.0;
    while ((j0 < x) | (j1 < x) | (j2 < n - x))
    {
        if ((j0 < x) && (f < 1))
        {
            j0++;
            f *= (double)(n - x + j0) / (double)j0;
        }
        else
        {
            if (j1 < x) 
            { 
                j1++; 
                f *= p; 
            }
            else 
            { 
                j2++; 
                f *= 1 - p; 
            }
        }
    }
    return(f);
}

#endif

#endif

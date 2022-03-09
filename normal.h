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
    return NAN; // x-mu is NaN 

  if (sigma == 0)
    return (x == mu) ? ML_POSINF : 0.;
    
  double x_ = (x - mu) / sigma;

  if (!isfinite(x_))
    return 0.;

  x_ = fabs(x_);
  if (x_ >= 2 * sqrt(DBL_MAX))
    return 0.;

  //  if (give_log)
  //    return -(M_LN_SQRT_2PI + 0.5 * x_ * x_ + log(sigma));
  //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)
//#ifdef MATHLIB_FAST_dnorm
  // and for R <= 3.0.x and R-devel upto 2014-01-01:
  return M_1_SQRT_2PI * exp(-0.5 * x_ * x_) / sigma;
//#else
/*
  // more accurate, less fast:
  if (x_ < 5)
    return M_1_SQRT_2PI * exp(-0.5 * x_ * x_) / sigma;
  // ELSE:
  // x*x  may lose upto about two digits accuracy for "large" x
  // Morten Welinder's proposal for PR#15620
  // https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620
  // -- 1 --  No hoop jumping when we underflow to zero anyway:
  //  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
  //     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
  // but "thanks" to denormalized numbers, underflow happens a bit later,
  //  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
  // for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
  // ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
  //              =IEEE=  38.58601
  // [on one x86_64 platform, effective boundary a bit lower: 38.56804]
  if (x_ > sqrt(-2 * M_LN2 * (DBL_MIN_EXP + 1 - DBL_MANT_DIG)))
    return 0.;
  // Now, to get full accurary, split x into two parts,
  //  x = x1+x2, such that |x2| <= 2^-16.
  // Assuming that we are using IEEE doubles, that means that
  // x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).
  // If we do not have IEEE this is still an improvement over the naive formula.
  double x1 = // R_forceint(x_ * 65536) / 65536 = ldexp(R_forceint(ldexp(x_, 16)), -16);
  double x2 = x_ - x1;
  return M_1_SQRT_2PI / sigma * (exp(-0.5 * x1 * x1) * exp((-0.5 * x2 - x1) * x2));
*/
//#endif
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

    // Handbook of Mathematical Functions by Abramowitz & Stegun, formula 7.1.26.
    double t = 1. / (1. + p * x);
    double y = 1. - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1. + sign * y);
}

/*
 *     Compute the quantile function for the normal distribution.
 *     For small to moderate probabilities, algorithm referenced
 *     below is used to obtain an initial approximation which is
 *     polished with a final Newton step.
 *     For very large arguments, an algorithm of Wichura is used.
 *  REFERENCE
 *     Beasley, J. D. and S. G. Springer (1977).
 *     Algorithm AS 111: The percentage points of the normal distribution,
 *     Applied Statistics, 26, 118-121.
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */
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

    ///-- use AS 241 --- 
    // double ppnd16_(double *p, long *ifault)
    //      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
    //      Produces the normal deviate Z corresponding to a given lower
    //      tail area of P; Z is accurate to about 1 part in 10**16.
    if (fabs(q) <= .425)
    {   // 0.075 <= p <= 0.925 
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
    {   // closer than 0.075 from {0,1} boundary 
        // r = min(p, 1-p) < 0.075 
        if (q > 0)
            r = 1 - p;
        else
            r = p;

        r = sqrt(-log(r));
        // r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) 

        if (r <= 5)
        {   // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 
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
        {   // very close to  0 or 1 
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

/*
// Implementation of the normal CDF with double precision approximations.
double phi(double x)
{
    static const double RT2PI = sqrt(4.0 * acos(0.0));

    static const double SPLIT = 7.07106781186547;

    static const double N0 = 220.206867912376;
    static const double N1 = 221.213596169931;
    static const double N2 = 112.079291497871;
    static const double N3 = 33.912866078383;
    static const double N4 = 6.37396220353165;
    static const double N5 = 0.700383064443688;
    static const double N6 = 3.52624965998911e-02;
    static const double M0 = 440.413735824752;
    static const double M1 = 793.826512519948;
    static const double M2 = 637.333633378831;
    static const double M3 = 296.564248779674;
    static const double M4 = 86.7807322029461;
    static const double M5 = 16.064177579207;
    static const double M6 = 1.75566716318264;
    static const double M7 = 8.83883476483184e-02;

    const double z = fabs(x);
    double c = 0.0;

    if (z <= 37.0)
    {
        const double e = exp(-z * z / 2.0);
        if (z < SPLIT)
        {
            const double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
            const double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
            c = e * n / d;
        }

        else
        {
            const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
            c = e / (RT2PI * f);
        }
    }

    return x <= 0.0 ? c : 1 - c;
}
*/

#endif

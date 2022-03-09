#ifndef CHISQUARE_H
#define CHISQUARE_H

#include "common.h"

// Computes the density of the gamma distribution,
double dgamma(double x, double shape, double scale)
{
    double pr;

#ifdef IEEE_754
    if (isnan(x) || isnan(shape) || isnan(scale))
        return x + shape + scale;
#endif

    if (shape < 0 || scale <= 0)
        return NAN;

    if (x < 0)
        return 0.;

    if (shape == 0) // point mass at 0 
        return (x == 0) ? ML_POSINF : 0.;

    if (x == 0)
    {
        if (shape < 1)
            return ML_POSINF;

        if (shape > 1)
            return 0.;

        // else 
        //return give_log ? -log(scale) : 1 / scale;
        return 1 / scale;
    }

    if (shape < 1)
    {
        pr = dpois_raw(shape, x / scale);
        // NB: currently *always*  shape/x > 0  if shape < 1:
        // -- overflow to Inf happens, but underflow to 0 does NOT: 
        // shape/x overflows to +Inf 
        //return (give_log ? pr + (isfinite(shape / x) ? log(shape / x) : log(shape) - log(x)) : pr * shape / x);
        return (pr * shape / x);
    }
    // else  shape >= 1 
    pr = dpois_raw(shape - 1, x / scale);
    //return give_log ? pr - log(scale) : pr / scale;
    return pr / scale;
}

// Compute the quantile function of the gamma distribution.
double qchisq_appr(double p, double nu, double g /* = log Gamma(nu/2) */, double tol /* EPS1 */)
{
#define C7	4.67
#define C8	6.66
#define C9	6.73
#define C10	13.32

    double alpha, a, c, ch, p1;
    double p2, q, t, x;

    // test arguments and initialise 

#ifdef IEEE_754
    if (isnan(p) || isnan(nu))
        return p + nu;
#endif

    if (p < 0 || p > 1)
        return NAN;

    if (nu <= 0)
        return NAN;

    alpha = 0.5 * nu;  // = [pq]gamma() shape
    c = alpha - 1;

    if (nu < (-1.24) * (p1 = log(p)))
    {	
        // for small chi-squared 
        // log(alpha) + g = log(alpha) + log(gamma(alpha)) =
        //   = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
        //   catastrophic cancellation when alpha << 1
        double lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g);
        ch = exp((lgam1pa + p1) / alpha + M_LN2);
    }
    else if (nu > 0.32)
    {
        //  using Wilson and Hilferty estimate 
        //x = qnorm(p, 0, 1, lower_tail, log_p);
        x = qNorm(p, 0., 1.);
        p1 = 2. / (9 * nu);
        ch = nu * pow(x * sqrt(p1) + 1 - p1, 3);

        // approximation for p tending to 1: 
        if (ch > 2.2 * nu + 6)
            ch = -2 * (log1p(-p) - c * log(0.5 * ch) + g);

    }
    else
    {
        // "small nu" : 1.24*(-log(p)) <= nu <= 0.32 
        ch = 0.4;
        a = log1p(p) + g + c * M_LN2;
        do {
            q = ch;
            p1 = 1. / (1 + ch * (C7 + ch));
            p2 = ch * (C9 + ch * (C8 + ch));
            t = -0.5 + (C7 + 2 * ch) * p1 - (C9 + ch * (C10 + 3 * ch)) / p2;
            ch -= (1 - exp(a + 0.5 * ch) * p2 * p1) / t;
        } while (fabs(q - ch) > tol * fabs(ch));
    }

    return ch;
}

double qgamma(double p, double alpha, double scale, int lower_tail) // shape = alpha 
{
#define EPS1   1e-2
#define EPS2   5e-7                // final precision of AS 91 
#define MAXIT  1000                // was 20 
#define pMIN   1e-100              // was 0.000002 = 2e-6 
#define pMAX   (1-1e-14)           // was (1-1e-12) and 0.999998 = 1 - 2e-6 

    const static double i420 = 1. / 420., i2520 = 1. / 2520., i5040 = 1. / 5040;
    double p_, a, b, c, g, ch, ch0, p1;
    double p2, q, s1, s2, s3, s4, s5, s6, t, x;

    // test arguments and initialise 

#ifdef IEEE_754
    if (isnan(p) || isnan(alpha) || isnan(scale))
        return p + alpha + scale;
#endif

    if (p < 0 || p > 1)
      return NAN;
    if (p == 0)
      return 0.;
    if (p == 1)
      return ML_POSINF;

    if (alpha < 0 || scale <= 0)
        return NAN;
    if (alpha == 0) // all mass at 0:
        return 0.;

    p_ = p; // lower_tail prob (in any case) 

    g = lgamma(alpha); // log Gamma(v/2) 

    //----- Phase I : Starting Approximation 

    ch = qchisq_appr(p, 2 * alpha, g, EPS1);
    if (!isfinite(ch))
        goto END;

    if (ch < EPS2)
        goto END;  // and do Newton steps 

    if (p_ > pMAX || p_ < pMIN)
        goto END;  // and do Newton steps 

      // ----- Phase II: Iteration
      // Call pgamma() [AS 239]	and calculate seven term taylor series

    c = alpha - 1;
    s6 = (120 + c * (346 + 127 * c)) * i5040;  // used below, is "const" 

    ch0 = ch;  // save initial approx. 
    for (unsigned i = 1; i <= MAXIT; i++)
    {
        q = ch;
        p1 = 0.5 * ch;
        p2 = p_ - pgamma_raw(p1, alpha, true);

#ifdef IEEE_754
        if (!isfinite(p2) || ch <= 0)
#else
        if (errno != 0 || ch <= 0)
#endif
        {
            ch = ch0;
            goto END;
        }
        
        t = p2 * exp(alpha * M_LN2 + g + p1 - c * log(ch));
        b = t / ch;
        a = 0.5 * t - b * c;
        s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) * i420;
        s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) * i2520;
        s3 = (210 + a * (462 + a * (707 + 932 * a))) * i2520;
        s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) * i5040;
        s5 = (84 + 2264 * a + c * (1175 + 606 * a)) * i2520;

        ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));

        if (fabs(q - ch) < EPS2 * ch)
            goto END;

        if (fabs(q - ch) > 0.1 * ch)
        {
            // diverging? -- also forces ch > 0 
            if (ch < q)
                ch = 0.9 * q; else ch = 1.1 * q;
        }
    }

END:
    x = 0.5 * scale * ch;

    return x;
}

// The density of the chi-squared distribution.
double dchisq(double x, double df) { return dgamma(x, df / 2., 2.); }

// The distribution function of the chi - squared distribution.
double pchisq(double x, double df) { return pgamma(x, df / 2., 2., true); }

// The quantile function of the chi-squared distribution.
double qchisq(double p, double df) { return qgamma(p, 0.5 * df, 2.0, true); }

#endif

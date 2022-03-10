#ifndef MACROS_H
#define MACROS_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <float.h>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cassert>

#define IEEE_754 1

extern double pNorm(const double, const double, const double);
extern double dNorm(const double x, const double, const double);
extern double qNorm(const double, const double, const double);

//#define ML_POSINF  std::numeric_limits<double>::infinity()
#define ML_POSINF  DBL_MAX
//#define ML_NEGINF  (-1 * std::numeric_limits<float>::infinity())
#define ML_NEGINF  (-1.0 * DBL_MAX)
#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI   0.398942280401432677939946059934 // 1/sqrt(2pi) 
#endif
#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif
#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI  0.918938533204672741780329736406 // log(sqrt(2*pi)) == log(2*pi)/2 
#endif
#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

// Utilities for use with the central limit theorem and normal distributions.
double qSigmaCLT(const double n, const double sigma) { return(sigma / sqrt(n)); }
double pSigmaCLT(const double n, const double p) { return (sqrt((p * (1. - p)) / n)); }
double zCLT(const double x, const double mu, const double sigma) { return ((x - mu) / sigma); }
double xCLT(const double z, const double mu, const double sigma) { return (mu + (z * sigma)); }
double zPhat(const double n, const double p, const double phat) { return (phat - p) / sqrt(p * (1. - p) / n); }

// Confidence intervals Margin of Error.
double proportionMoE(const double n, const double z, const double phat) { return (z * sqrt(phat * (1. - phat) / n)); }
double meanMoE(const double n, const double t, const double sigma) { return (t * (sigma / sqrt(n))); }
// Confidence intervals n. 
double proportionN(const double MoE, const double z, const double phat) { return (phat * (1. - phat) * pow((z / MoE), 2.)); }
double meanN(const double MoE, const double z, const double sigma) { return pow((z * sigma) / MoE, 2); }
// Confidence intervals z-scores.
constexpr double Z95CI = 1.95996; // -1.95996 (right tail). Z95CI = qNorm(.95 + (1 - .95)/ 2);
constexpr double Z90CI = 1.64485; // -1.64485 (right tail).

// Hypothesis testing z and t.
double proportionHypothesisZ(const unsigned n, const double phat, const double p0) { return ((phat - p0) / sqrt((p0 * (1. - p0)) / n)); }
double meanHypothesisT(const unsigned n, const double xbar, const double mu, const double sigma) { return ((xbar - mu) / (sigma / sqrt(n))); }
#define DecideHypothesis(pValue, alpha)  (pValue < alpha) ? std::cout << "reject H0\n" : std::cout << "don't reject H0\n"

double fmax2(double x, double y)
{
#ifdef IEEE_754
    if (isnan(x) || isnan(y))
        return x + y;
#endif
    return (x < y) ? y : x;
}

double fmin2(double x, double y)
{
#ifdef IEEE_754
    if (isnan(x) || isnan(y))
        return x + y;
#endif
    return (x < y) ? x : y;
}

// Computes the log of the error term in Stirling's formula.
double stirlerr(double n)
{
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
    // error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
    const static double sferr_halves[31] = {
    0.0, /* n=0 - wrong, place holder only */
    0.1534264097200273452913848,  /* 0.5 */
    0.0810614667953272582196702,  /* 1.0 */
    0.0548141210519176538961390,  /* 1.5 */
    0.0413406959554092940938221,  /* 2.0 */
    0.03316287351993628748511048, /* 2.5 */
    0.02767792568499833914878929, /* 3.0 */
    0.02374616365629749597132920, /* 3.5 */
    0.02079067210376509311152277, /* 4.0 */
    0.01848845053267318523077934, /* 4.5 */
    0.01664469118982119216319487, /* 5.0 */
    0.01513497322191737887351255, /* 5.5 */
    0.01387612882307074799874573, /* 6.0 */
    0.01281046524292022692424986, /* 6.5 */
    0.01189670994589177009505572, /* 7.0 */
    0.01110455975820691732662991, /* 7.5 */
    0.010411265261972096497478567, /* 8.0 */
    0.009799416126158803298389475, /* 8.5 */
    0.009255462182712732917728637, /* 9.0 */
    0.008768700134139385462952823, /* 9.5 */
    0.008330563433362871256469318, /* 10.0 */
    0.007934114564314020547248100, /* 10.5 */
    0.007573675487951840794972024, /* 11.0 */
    0.007244554301320383179543912, /* 11.5 */
    0.006942840107209529865664152, /* 12.0 */
    0.006665247032707682442354394, /* 12.5 */
    0.006408994188004207068439631, /* 13.0 */
    0.006171712263039457647532867, /* 13.5 */
    0.005951370112758847735624416, /* 14.0 */
    0.005746216513010115682023589, /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0)
    {
        nn = n + n;
        if (nn == (int)nn)
            return(sferr_halves[(int)nn]);
        //return(lgammafn(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI);
        return(lgamma(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n * n;

    if (n > 500)
        return((S0 - S1 / nn) / n);

    if (n > 80)
        return((S0 - (S1 - S2 / nn) / nn) / n);

    if (n > 35)
        return((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);

    // 15 < n <= 35 : 
    return((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
}

#define M_SQRT_2PI  2.50662827463100050241576528481104525301  // sqrt(2*pi) 
// sqrt(2 * Rmpfr::Const("pi", 128))
#define x_LRG  2.86111748575702815380240589208115399625e+307  // = 2^1023 / pi
// Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

static double logcf(double x, double i, double d, double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

    assert(i > 0);
    assert(d >= 0);

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2))
    {
        double c3 = c2 * c2 * x;

        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;

        c3 = c1 * c1 * x;
        c1 += d;
        c4 += d;
        a2 = c4 * a1 - c3 * a2;
        b2 = c4 * b1 - c3 * b2;

        if (fabs(b2) > scalefactor)
        {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        }
        else if (fabs(b2) < 1 / scalefactor)
        {
            a1 *= scalefactor;
            b1 *= scalefactor;
            a2 *= scalefactor;
            b2 *= scalefactor;
        }
    }

    return a2 / b2;
}

// Accurate calculation of log(1+x)-x, particularly for small x.  
double log1pmx(double x) {
    static const double minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
        return log1p(x) - x;
    else
    {
        double r = x / (2 + x), y = r * r;

        if (fabs(x) < 1e-2)
        {
            static const double two = 2;
            return r * ((((two / 9 * y + two / 7) * y + two / 5) * y + two / 3) * y - x);
        }
        else
        {
            static const double tol_logcf = 1e-14;
            return r * (2 * y * logcf(y, 3, 2, tol_logcf) - x);
        }
    }
}

static const float bd0_scale[128 + 1][4] = {
  { +0x1.62e430p-1, -0x1.05c610p-29, -0x1.950d88p-54, +0x1.d9cc02p-79 }, // 128: log(2048/1024.) 
  { +0x1.5ee02cp-1, -0x1.6dbe98p-25, -0x1.51e540p-50, +0x1.2bfa48p-74 }, // 129: log(2032/1024.) 
  { +0x1.5ad404p-1, +0x1.86b3e4p-26, +0x1.9f6534p-50, +0x1.54be04p-74 }, // 130: log(2016/1024.) 
  { +0x1.570124p-1, -0x1.9ed750p-25, -0x1.f37dd0p-51, +0x1.10b770p-77 }, // 131: log(2001/1024.) 
  { +0x1.5326e4p-1, -0x1.9b9874p-25, -0x1.378194p-49, +0x1.56feb2p-74 }, // 132: log(1986/1024.) 
  { +0x1.4f4528p-1, +0x1.aca70cp-28, +0x1.103e74p-53, +0x1.9c410ap-81 }, // 133: log(1971/1024.) 
  { +0x1.4b5bd8p-1, -0x1.6a91d8p-25, -0x1.8e43d0p-50, -0x1.afba9ep-77 }, // 134: log(1956/1024.) 
  { +0x1.47ae54p-1, -0x1.abb51cp-25, +0x1.19b798p-51, +0x1.45e09cp-76 }, // 135: log(1942/1024.) 
  { +0x1.43fa00p-1, -0x1.d06318p-25, -0x1.8858d8p-49, -0x1.1927c4p-75 }, // 136: log(1928/1024.) 
  { +0x1.3ffa40p-1, +0x1.1a427cp-25, +0x1.151640p-53, -0x1.4f5606p-77 }, // 137: log(1913/1024.) 
  { +0x1.3c7c80p-1, -0x1.19bf48p-34, +0x1.05fc94p-58, -0x1.c096fcp-82 }, // 138: log(1900/1024.) 
  { +0x1.38b320p-1, +0x1.6b5778p-25, +0x1.be38d0p-50, -0x1.075e96p-74 }, // 139: log(1886/1024.) 
  { +0x1.34e288p-1, +0x1.d9ce1cp-25, +0x1.316eb8p-49, +0x1.2d885cp-73 }, // 140: log(1872/1024.) 
  { +0x1.315124p-1, +0x1.c2fc60p-29, -0x1.4396fcp-53, +0x1.acf376p-78 }, // 141: log(1859/1024.) 
  { +0x1.2db954p-1, +0x1.720de4p-25, -0x1.d39b04p-49, -0x1.f11176p-76 }, // 142: log(1846/1024.) 
  { +0x1.2a1b08p-1, -0x1.562494p-25, +0x1.a7863cp-49, +0x1.85dd64p-73 }, // 143: log(1833/1024.) 
  { +0x1.267620p-1, +0x1.3430e0p-29, -0x1.96a958p-56, +0x1.f8e636p-82 }, // 144: log(1820/1024.) 
  { +0x1.23130cp-1, +0x1.7bebf4p-25, +0x1.416f1cp-52, -0x1.78dd36p-77 }, // 145: log(1808/1024.) 
  { +0x1.1faa34p-1, +0x1.70e128p-26, +0x1.81817cp-50, -0x1.c2179cp-76 }, // 146: log(1796/1024.) 
  { +0x1.1bf204p-1, +0x1.3a9620p-28, +0x1.2f94c0p-52, +0x1.9096c0p-76 }, // 147: log(1783/1024.) 
  { +0x1.187ce4p-1, -0x1.077870p-27, +0x1.655a80p-51, +0x1.eaafd6p-78 }, // 148: log(1771/1024.) 
  { +0x1.1501c0p-1, -0x1.406cacp-25, -0x1.e72290p-49, +0x1.5dd800p-73 }, // 149: log(1759/1024.) 
  { +0x1.11cb80p-1, +0x1.787cd0p-25, -0x1.efdc78p-51, -0x1.5380cep-77 }, // 150: log(1748/1024.) 
  { +0x1.0e4498p-1, +0x1.747324p-27, -0x1.024548p-51, +0x1.77a5a6p-75 }, // 151: log(1736/1024.) 
  { +0x1.0b036cp-1, +0x1.690c74p-25, +0x1.5d0cc4p-50, -0x1.c0e23cp-76 }, // 152: log(1725/1024.) 
  { +0x1.077070p-1, -0x1.a769bcp-27, +0x1.452234p-52, +0x1.6ba668p-76 }, // 153: log(1713/1024.) 
  { +0x1.04240cp-1, -0x1.a686acp-27, -0x1.ef46b0p-52, -0x1.5ce10cp-76 }, // 154: log(1702/1024.) 
  { +0x1.00d22cp-1, +0x1.fc0e10p-25, +0x1.6ee034p-50, -0x1.19a2ccp-74 }, // 155: log(1691/1024.) 
  { +0x1.faf588p-2, +0x1.ef1e64p-27, -0x1.26504cp-54, -0x1.b15792p-82 }, // 156: log(1680/1024.) 
  { +0x1.f4d87cp-2, +0x1.d7b980p-26, -0x1.a114d8p-50, +0x1.9758c6p-75 }, // 157: log(1670/1024.) 
  { +0x1.ee1414p-2, +0x1.2ec060p-26, +0x1.dc00fcp-52, +0x1.f8833cp-76 }, // 158: log(1659/1024.) 
  { +0x1.e7e32cp-2, -0x1.ac796cp-27, -0x1.a68818p-54, +0x1.235d02p-78 }, // 159: log(1649/1024.) 
  { +0x1.e108a0p-2, -0x1.768ba4p-28, -0x1.f050a8p-52, +0x1.00d632p-82 }, // 160: log(1638/1024.) 
  { +0x1.dac354p-2, -0x1.d3a6acp-30, +0x1.18734cp-57, -0x1.f97902p-83 }, // 161: log(1628/1024.) 
  { +0x1.d47424p-2, +0x1.7dbbacp-31, -0x1.d5ada4p-56, +0x1.56fcaap-81 }, // 162: log(1618/1024.) 
  { +0x1.ce1af0p-2, +0x1.70be7cp-27, +0x1.6f6fa4p-51, +0x1.7955a2p-75 }, // 163: log(1608/1024.) 
  { +0x1.c7b798p-2, +0x1.ec36ecp-26, -0x1.07e294p-50, -0x1.ca183cp-75 }, // 164: log(1598/1024.) 
  { +0x1.c1ef04p-2, +0x1.c1dfd4p-26, +0x1.888eecp-50, -0x1.fd6b86p-75 }, // 165: log(1589/1024.) 
  { +0x1.bb7810p-2, +0x1.478bfcp-26, +0x1.245b8cp-50, +0x1.ea9d52p-74 }, // 166: log(1579/1024.) 
  { +0x1.b59da0p-2, -0x1.882b08p-27, +0x1.31573cp-53, -0x1.8c249ap-77 }, // 167: log(1570/1024.) 
  { +0x1.af1294p-2, -0x1.b710f4p-27, +0x1.622670p-51, +0x1.128578p-76 }, // 168: log(1560/1024.) 
  { +0x1.a925d4p-2, -0x1.0ae750p-27, +0x1.574ed4p-51, +0x1.084996p-75 }, // 169: log(1551/1024.) 
  { +0x1.a33040p-2, +0x1.027d30p-29, +0x1.b9a550p-53, -0x1.b2e38ap-78 }, // 170: log(1542/1024.) 
  { +0x1.9d31c0p-2, -0x1.5ec12cp-26, -0x1.5245e0p-52, +0x1.2522d0p-79 }, // 171: log(1533/1024.) 
  { +0x1.972a34p-2, +0x1.135158p-30, +0x1.a5c09cp-56, +0x1.24b70ep-80 }, // 172: log(1524/1024.) 
  { +0x1.911984p-2, +0x1.0995d4p-26, +0x1.3bfb5cp-50, +0x1.2c9dd6p-75 }, // 173: log(1515/1024.) 
  { +0x1.8bad98p-2, -0x1.1d6144p-29, +0x1.5b9208p-53, +0x1.1ec158p-77 }, // 174: log(1507/1024.) 
  { +0x1.858b58p-2, -0x1.1b4678p-27, +0x1.56cab4p-53, -0x1.2fdc0cp-78 }, // 175: log(1498/1024.) 
  { +0x1.7f5fa0p-2, +0x1.3aaf48p-27, +0x1.461964p-51, +0x1.4ae476p-75 }, // 176: log(1489/1024.) 
  { +0x1.79db68p-2, -0x1.7e5054p-26, +0x1.673750p-51, -0x1.a11f7ap-76 }, // 177: log(1481/1024.) 
  { +0x1.744f88p-2, -0x1.cc0e18p-26, -0x1.1e9d18p-50, -0x1.6c06bcp-78 }, // 178: log(1473/1024.) 
  { +0x1.6e08ecp-2, -0x1.5d45e0p-26, -0x1.c73ec8p-50, +0x1.318d72p-74 }, // 179: log(1464/1024.) 
  { +0x1.686c80p-2, +0x1.e9b14cp-26, -0x1.13bbd4p-50, -0x1.efeb1cp-78 }, // 180: log(1456/1024.) 
  { +0x1.62c830p-2, -0x1.a8c70cp-27, -0x1.5a1214p-51, -0x1.bab3fcp-79 }, // 181: log(1448/1024.) 
  { +0x1.5d1bdcp-2, -0x1.4fec6cp-31, +0x1.423638p-56, +0x1.ee3feep-83 }, // 182: log(1440/1024.) 
  { +0x1.576770p-2, +0x1.7455a8p-26, -0x1.3ab654p-50, -0x1.26be4cp-75 }, // 183: log(1432/1024.) 
  { +0x1.5262e0p-2, -0x1.146778p-26, -0x1.b9f708p-52, -0x1.294018p-77 }, // 184: log(1425/1024.) 
  { +0x1.4c9f08p-2, +0x1.e152c4p-26, -0x1.dde710p-53, +0x1.fd2208p-77 }, // 185: log(1417/1024.) 
  { +0x1.46d2d8p-2, +0x1.c28058p-26, -0x1.936284p-50, +0x1.9fdd68p-74 }, // 186: log(1409/1024.) 
  { +0x1.41b940p-2, +0x1.cce0c0p-26, -0x1.1a4050p-50, +0x1.bc0376p-76 }, // 187: log(1402/1024.) 
  { +0x1.3bdd24p-2, +0x1.d6296cp-27, +0x1.425b48p-51, -0x1.cddb2cp-77 }, // 188: log(1394/1024.) 
  { +0x1.36b578p-2, -0x1.287ddcp-27, -0x1.2d0f4cp-51, +0x1.38447ep-75 }, // 189: log(1387/1024.) 
  { +0x1.31871cp-2, +0x1.2a8830p-27, +0x1.3eae54p-52, -0x1.898136p-77 }, // 190: log(1380/1024.) 
  { +0x1.2b9304p-2, -0x1.51d8b8p-28, +0x1.27694cp-52, -0x1.fd852ap-76 }, // 191: log(1372/1024.) 
  { +0x1.265620p-2, -0x1.d98f3cp-27, +0x1.a44338p-51, -0x1.56e85ep-78 }, // 192: log(1365/1024.) 
  { +0x1.211254p-2, +0x1.986160p-26, +0x1.73c5d0p-51, +0x1.4a861ep-75 }, // 193: log(1358/1024.) 
  { +0x1.1bc794p-2, +0x1.fa3918p-27, +0x1.879c5cp-51, +0x1.16107cp-78 }, // 194: log(1351/1024.) 
  { +0x1.1675ccp-2, -0x1.4545a0p-26, +0x1.c07398p-51, +0x1.f55c42p-76 }, // 195: log(1344/1024.) 
  { +0x1.111ce4p-2, +0x1.f72670p-37, -0x1.b84b5cp-61, +0x1.a4a4dcp-85 }, // 196: log(1337/1024.) 
  { +0x1.0c81d4p-2, +0x1.0c150cp-27, +0x1.218600p-51, -0x1.d17312p-76 }, // 197: log(1331/1024.) 
  { +0x1.071b84p-2, +0x1.fcd590p-26, +0x1.a3a2e0p-51, +0x1.fe5ef8p-76 }, // 198: log(1324/1024.) 
  { +0x1.01ade4p-2, -0x1.bb1844p-28, +0x1.db3cccp-52, +0x1.1f56fcp-77 }, // 199: log(1317/1024.) 
  { +0x1.fa01c4p-3, -0x1.12a0d0p-29, -0x1.f71fb0p-54, +0x1.e287a4p-78 }, // 200: log(1311/1024.) 
  { +0x1.ef0adcp-3, +0x1.7b8b28p-28, -0x1.35bce4p-52, -0x1.abc8f8p-79 }, // 201: log(1304/1024.) 
  { +0x1.e598ecp-3, +0x1.5a87e4p-27, -0x1.134bd0p-51, +0x1.c2cebep-76 }, // 202: log(1298/1024.) 
  { +0x1.da85d8p-3, -0x1.df31b0p-27, +0x1.94c16cp-57, +0x1.8fd7eap-82 }, // 203: log(1291/1024.) 
  { +0x1.d0fb80p-3, -0x1.bb5434p-28, -0x1.ea5640p-52, -0x1.8ceca4p-77 }, // 204: log(1285/1024.) 
  { +0x1.c765b8p-3, +0x1.e4d68cp-27, +0x1.5b59b4p-51, +0x1.76f6c4p-76 }, // 205: log(1279/1024.) 
  { +0x1.bdc46cp-3, -0x1.1cbb50p-27, +0x1.2da010p-51, +0x1.eb282cp-75 }, // 206: log(1273/1024.) 
  { +0x1.b27980p-3, -0x1.1b9ce0p-27, +0x1.7756f8p-52, +0x1.2ff572p-76 }, // 207: log(1266/1024.) 
  { +0x1.a8bed0p-3, -0x1.bbe874p-30, +0x1.85cf20p-56, +0x1.b9cf18p-80 }, // 208: log(1260/1024.) 
  { +0x1.9ef83cp-3, +0x1.2769a4p-27, -0x1.85bda0p-52, +0x1.8c8018p-79 }, // 209: log(1254/1024.) 
  { +0x1.9525a8p-3, +0x1.cf456cp-27, -0x1.7137d8p-52, -0x1.f158e8p-76 }, // 210: log(1248/1024.) 
  { +0x1.8b46f8p-3, +0x1.11b12cp-30, +0x1.9f2104p-54, -0x1.22836ep-78 }, // 211: log(1242/1024.) 
  { +0x1.83040cp-3, +0x1.2379e4p-28, +0x1.b71c70p-52, -0x1.990cdep-76 }, // 212: log(1237/1024.) 
  { +0x1.790ed4p-3, +0x1.dc4c68p-28, -0x1.910ac8p-52, +0x1.dd1bd6p-76 }, // 213: log(1231/1024.) 
  { +0x1.6f0d28p-3, +0x1.5cad68p-28, +0x1.737c94p-52, -0x1.9184bap-77 }, // 214: log(1225/1024.) 
  { +0x1.64fee8p-3, +0x1.04bf88p-28, +0x1.6fca28p-52, +0x1.8884a8p-76 }, // 215: log(1219/1024.) 
  { +0x1.5c9400p-3, +0x1.d65cb0p-29, -0x1.b2919cp-53, +0x1.b99bcep-77 }, // 216: log(1214/1024.) 
  { +0x1.526e60p-3, -0x1.c5e4bcp-27, -0x1.0ba380p-52, +0x1.d6e3ccp-79 }, // 217: log(1208/1024.) 
  { +0x1.483bccp-3, +0x1.9cdc7cp-28, -0x1.5ad8dcp-54, -0x1.392d3cp-83 }, // 218: log(1202/1024.) 
  { +0x1.3fb25cp-3, -0x1.a6ad74p-27, +0x1.5be6b4p-52, -0x1.4e0114p-77 }, // 219: log(1197/1024.) 
  { +0x1.371fc4p-3, -0x1.fe1708p-27, -0x1.78864cp-52, -0x1.27543ap-76 }, // 220: log(1192/1024.) 
  { +0x1.2cca10p-3, -0x1.4141b4p-28, -0x1.ef191cp-52, +0x1.00ee08p-76 }, // 221: log(1186/1024.) 
  { +0x1.242310p-3, +0x1.3ba510p-27, -0x1.d003c8p-51, +0x1.162640p-76 }, // 222: log(1181/1024.) 
  { +0x1.1b72acp-3, +0x1.52f67cp-27, -0x1.fd6fa0p-51, +0x1.1a3966p-77 }, // 223: log(1176/1024.) 
  { +0x1.10f8e4p-3, +0x1.129cd8p-30, +0x1.31ef30p-55, +0x1.a73e38p-79 }, // 224: log(1170/1024.) 
  { +0x1.08338cp-3, -0x1.005d7cp-27, -0x1.661a9cp-51, +0x1.1f138ap-79 }, // 225: log(1165/1024.) 
  { +0x1.fec914p-4, -0x1.c482a8p-29, -0x1.55746cp-54, +0x1.99f932p-80 }, // 226: log(1160/1024.) 
  { +0x1.ed1794p-4, +0x1.d06f00p-29, +0x1.75e45cp-53, -0x1.d0483ep-78 }, // 227: log(1155/1024.) 
  { +0x1.db5270p-4, +0x1.87d928p-32, -0x1.0f52a4p-57, +0x1.81f4a6p-84 }, // 228: log(1150/1024.) 
  { +0x1.c97978p-4, +0x1.af1d24p-29, -0x1.0977d0p-60, -0x1.8839d0p-84 }, // 229: log(1145/1024.) 
  { +0x1.b78c84p-4, -0x1.44f124p-28, -0x1.ef7bc4p-52, +0x1.9e0650p-78 }, // 230: log(1140/1024.) 
  { +0x1.a58b60p-4, +0x1.856464p-29, +0x1.c651d0p-55, +0x1.b06b0cp-79 }, // 231: log(1135/1024.) 
  { +0x1.9375e4p-4, +0x1.5595ecp-28, +0x1.dc3738p-52, +0x1.86c89ap-81 }, // 232: log(1130/1024.) 
  { +0x1.814be4p-4, -0x1.c073fcp-28, -0x1.371f88p-53, -0x1.5f4080p-77 }, // 233: log(1125/1024.) 
  { +0x1.6f0d28p-4, +0x1.5cad68p-29, +0x1.737c94p-53, -0x1.9184bap-78 }, // 234: log(1120/1024.) 
  { +0x1.60658cp-4, -0x1.6c8af4p-28, +0x1.d8ef74p-55, +0x1.c4f792p-80 }, // 235: log(1116/1024.) 
  { +0x1.4e0110p-4, +0x1.146b5cp-29, +0x1.73f7ccp-54, -0x1.d28db8p-79 }, // 236: log(1111/1024.) 
  { +0x1.3b8758p-4, +0x1.8b1b70p-28, -0x1.20aca4p-52, -0x1.651894p-76 }, // 237: log(1106/1024.) 
  { +0x1.28f834p-4, +0x1.43b6a4p-30, -0x1.452af8p-55, +0x1.976892p-80 }, // 238: log(1101/1024.) 
  { +0x1.1a0fbcp-4, -0x1.e4075cp-28, +0x1.1fe618p-52, +0x1.9d6dc2p-77 }, // 239: log(1097/1024.) 
  { +0x1.075984p-4, -0x1.4ce370p-29, -0x1.d9fc98p-53, +0x1.4ccf12p-77 }, // 240: log(1092/1024.) 
  { +0x1.f0a30cp-5, +0x1.162a68p-37, -0x1.e83368p-61, -0x1.d222a6p-86 }, // 241: log(1088/1024.) 
  { +0x1.cae730p-5, -0x1.1a8f7cp-31, -0x1.5f9014p-55, +0x1.2720c0p-79 }, // 242: log(1083/1024.) 
  { +0x1.ac9724p-5, -0x1.e8ee08p-29, +0x1.a7de04p-54, -0x1.9bba74p-78 }, // 243: log(1079/1024.) 
  { +0x1.868a84p-5, -0x1.ef8128p-30, +0x1.dc5eccp-54, -0x1.58d250p-79 }, // 244: log(1074/1024.) 
  { +0x1.67f950p-5, -0x1.ed684cp-30, -0x1.f060c0p-55, -0x1.b1294cp-80 }, // 245: log(1070/1024.) 
  { +0x1.494accp-5, +0x1.a6c890p-32, -0x1.c3ad48p-56, -0x1.6dc66cp-84 }, // 246: log(1066/1024.) 
  { +0x1.22c71cp-5, -0x1.8abe2cp-32, -0x1.7e7078p-56, -0x1.ddc3dcp-86 }, // 247: log(1061/1024.) 
  { +0x1.03d5d8p-5, +0x1.79cfbcp-31, -0x1.da7c4cp-58, +0x1.4e7582p-83 }, // 248: log(1057/1024.) 
  { +0x1.c98d18p-6, +0x1.a01904p-31, -0x1.854164p-55, +0x1.883c36p-79 }, // 249: log(1053/1024.) 
  { +0x1.8b31fcp-6, -0x1.356500p-30, +0x1.c3ab48p-55, +0x1.b69bdap-80 }, // 250: log(1049/1024.) 
  { +0x1.3cea44p-6, +0x1.a352bcp-33, -0x1.8865acp-57, -0x1.48159cp-81 }, // 251: log(1044/1024.) 
  { +0x1.fc0a8cp-7, -0x1.e07f84p-32, +0x1.e7cf6cp-58, +0x1.3a69c0p-82 }, // 252: log(1040/1024.) 
  { +0x1.7dc474p-7, +0x1.f810a8p-31, -0x1.245b5cp-56, -0x1.a1f4f8p-80 }, // 253: log(1036/1024.) 
  { +0x1.fe02a8p-8, -0x1.4ef988p-32, +0x1.1f86ecp-57, +0x1.20723cp-81 }, // 254: log(1032/1024.) 
  { +0x1.ff00acp-9, -0x1.d4ef44p-33, +0x1.2821acp-63, +0x1.5a6d32p-87 }, // 255: log(1028/1024.) 
  { 0, 0, 0, 0 }
};

#define ADD1(d_) do {              \
      double d = (d_);             \
      double d1 = floor (d + 0.5); \
      double d2 = d - d1;          \
      *yh += d1;                   \
      *yl += d2;                   \
  } while(0)

// Compute x * log (x / M) + (M - x)
void ebd0(double x, double M, double* yh, double* yl)
{
    const int Sb = 10;
    const double S = 1u << Sb;
    const int N = 128;

    *yl = *yh = 0;

    if (x == M)
        return;

    if (x == 0)
    {
        *yh = M;
        return;
    }

    if (M == 0)
    {
        *yh = ML_POSINF;
        return;
    }

    if (M / x == ML_POSINF)
    {
        *yh = M;
        return;
    }

    int e;
    double r = frexp(M / x, &e);

    // prevent later overflow
    if (M_LN2 * ((double)-e) > 1. + DBL_MAX / x)
    {
        *yh = ML_POSINF;
        return;
    }

    int i = (int)floor((r - 0.5) * (2 * N) + 0.5);
    // now,  0 <= i <= N
    double f = floor(S / (0.5 + i / (2.0 * N)) + 0.5);
    double fg = ldexp(f, -(e + Sb)); // ldexp(f, E) := f * 2^E

    if (fg == ML_POSINF)
    {
        *yh = fg;
        return;
    }

    ADD1(-x * log1pmx((M * fg - x) / x));
    if (fg == 1)
        return;
    // else (fg != 1) :
    for (int j = 0; j < 4; j++)
    {
        ADD1(x * bd0_scale[i][j]);  // handles  x*log(fg*2^e) 
        ADD1(-x * e * bd0_scale[0][j]);  // handles  x*log(1/ 2^e) 
        if (!isfinite(*yh))
        {
            *yh = ML_POSINF;
            *yl = 0; return;
        }
    }

    ADD1(M);
    ADD1(-M * fg);
}

#undef ADD1

//  dpois_raw() computes the Poisson probability  lb^x exp(-lb) / x!.
double dpois_raw(double x, double lambda)
{
    if (lambda == 0)
        return ((x == 0) ? 1. : 0.);

    if (!isfinite(lambda))
        return 0.; // including for the case where  x = lambda = +Inf

    if (x < 0.)
        return 0.;

    if (x <= lambda * DBL_MIN)
        return (exp(-lambda));

    if (lambda < x * DBL_MIN)
    {
        if (!isfinite(x)) // lambda < x = +Inf
            return 0.;
        // else
        return(exp(-lambda + x * log(lambda) - lgamma(x + 1)));
    }

    double yh, yl;

    ebd0(x, lambda, &yh, &yl);
    yl += stirlerr(x);

    bool Lrg_x = (x >= x_LRG); //really large x  <==>  2*pi*x  overflows
    double r = Lrg_x ? M_SQRT_2PI * sqrt(x) : M_2PI * x; // sqrt(.): avoid overflow for very large x

    return exp(-yl) * exp(-yh) / (Lrg_x ? r : sqrt(r));
}

// Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). 
double lgamma1p(double a)
{
    if (fabs(a) >= 0.5)
        //return lgammafn(a + 1);
        return lgamma(a + 1);

    const double eulers_const = 0.5772156649015328606065120900824024;

    // coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40:
    const int N = 40;
    static const double coeffs[40] = {
      0.3224670334241132182362075833230126e-0,  // = (zeta(2)-1)/2 
      0.6735230105319809513324605383715000e-1,  // = (zeta(3)-1)/3 
      0.2058080842778454787900092413529198e-1,
      0.7385551028673985266273097291406834e-2,
      0.2890510330741523285752988298486755e-2,
      0.1192753911703260977113935692828109e-2,
      0.5096695247430424223356548135815582e-3,
      0.2231547584535793797614188036013401e-3,
      0.9945751278180853371459589003190170e-4,
      0.4492623673813314170020750240635786e-4,
      0.2050721277567069155316650397830591e-4,
      0.9439488275268395903987425104415055e-5,
      0.4374866789907487804181793223952411e-5,
      0.2039215753801366236781900709670839e-5,
      0.9551412130407419832857179772951265e-6,
      0.4492469198764566043294290331193655e-6,
      0.2120718480555466586923135901077628e-6,
      0.1004322482396809960872083050053344e-6,
      0.4769810169363980565760193417246730e-7,
      0.2271109460894316491031998116062124e-7,
      0.1083865921489695409107491757968159e-7,
      0.5183475041970046655121248647057669e-8,
      0.2483674543802478317185008663991718e-8,
      0.1192140140586091207442548202774640e-8,
      0.5731367241678862013330194857961011e-9,
      0.2759522885124233145178149692816341e-9,
      0.1330476437424448948149715720858008e-9,
      0.6422964563838100022082448087644648e-10,
      0.3104424774732227276239215783404066e-10,
      0.1502138408075414217093301048780668e-10,
      0.7275974480239079662504549924814047e-11,
      0.3527742476575915083615072228655483e-11,
      0.1711991790559617908601084114443031e-11,
      0.8315385841420284819798357793954418e-12,
      0.4042200525289440065536008957032895e-12,
      0.1966475631096616490411045679010286e-12,
      0.9573630387838555763782200936508615e-13,
      0.4664076026428374224576492565974577e-13,
      0.2273736960065972320633279596737272e-13,
      0.1109139947083452201658320007192334e-13 // = (zeta(40+1)-1)/(40+1) 
    };

    const double c = 0.2273736845824652515226821577978691e-12; // zeta(N+2)-1 
    const double tol_logcf = 1e-14;

    double lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);

    for (int i = N - 1; i >= 0; i--)
        lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx(a);
}

// Abramowitz and Stegun 6.5.29 [right]
static double pgamma_smallx(double x, double alph, int lower_tail)
{
    double sum = 0, c = alph, n = 0, term;

    do {
        n++;
        c *= -x / n;
        term = c / (alph + n);
        sum += term;
    } while (fabs(term) > DBL_EPSILON * fabs(sum));

    if (lower_tail)
    {
        double f1 = 1 + sum;
        double f2;
        if (alph > 1)
        {
            f2 = dpois_raw(alph, x);
            f2 = f2 * exp(x);
        }
        else
            f2 = pow(x, alph) / exp(lgamma1p(alph));
        return f1 * f2;
    }
    else
    {
        double lf2 = alph * log(x) - lgamma1p(alph);
        double f1m1 = sum;
        double f2m1 = expm1(lf2);
        return -(f1m1 + f2m1 + f1m1 * f2m1);
    }
}

static double pd_upper_series(double x, double y)
{
    double term = x / y;
    double sum = term;

    do {
        y++;
        term *= x / y;
        sum += term;
    } while (term > sum * DBL_EPSILON);

    return sum;
}

// Continued fraction for calculation of scaled upper-tail F_{gamma}
static double pd_lower_cf(double y, double d)
{
    double f = 0.0 /* -Wall */, of, f0;
    double i, c2, c3, c4, a1, b1, a2, b2;

#define  NEEDED_SCALE   \
  (b2 > scalefactor) {  \
     a1 /= scalefactor; \
     b1 /= scalefactor; \
     a2 /= scalefactor; \
     b2 /= scalefactor; \
  }

#define max_it 200000

    if (y == 0)
        return 0;

    f0 = y / d;
    if (fabs(y - 1) < fabs(d) * DBL_EPSILON)
        return (f0);

    if (f0 > 1.)
        f0 = 1.;
    c2 = y;
    c4 = d;

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

        i = 0; of = -1.; // far away 
    while (i < max_it)
    {
        i++;  c2--;  c3 = i * c2;  c4 += 2;
        // c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd 
        a1 = c4 * a2 + c3 * a1;
        b1 = c4 * b2 + c3 * b1;

        i++;  c2--;  c3 = i * c2;  c4 += 2;
        // c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even 
        a2 = c4 * a1 + c3 * a2;
        b2 = c4 * b1 + c3 * b2;

        if NEEDED_SCALE

            if (b2 != 0)
            {
                f = a2 / b2;
                if (fabs(f - of) <= DBL_EPSILON * fmax2(f0, fabs(f)))
                    return f;
                of = f;
            }
    }

    return f; // should not happen... 
}

#undef NEEDED_SCALE

static double pd_lower_series(double lambda, double y)
{
    double term = 1, sum = 0;

    while (y >= 1 && term > sum * DBL_EPSILON)
    {
        term *= y / lambda;
        sum += term;
        y--;
    }

    if (y != floor(y))
    {
        double f;

        f = pd_lower_cf(y, lambda + 1 - y);
        sum += term * f;
    }

    return sum;
}

// Asymptotic expansion to calculate the probability that Poisson variate has value <= x.
static double ppois_asymp(double x, double lambda, bool lower_tail)
{
    static const double coefs_a[8] =
    {
      -1e99, // placeholder used for 1-indexing 
      2 / 3.,
      -4 / 135.,
      8 / 2835.,
      16 / 8505.,
      -8992 / 12629925.,
      -334144 / 492567075.,
      698752 / 1477701225.
    };

    static const double coefs_b[8] =
    {
      -1e99, // placeholder 
      1 / 12.,
      1 / 288.,
      -139 / 51840.,
      -571 / 2488320.,
      163879 / 209018880.,
      5246819 / 75246796800.,
      -534703531 / 902961561600.
    };

    double elfb, elfb_term;
    double res12, res1_term, res1_ig, res2_term, res2_ig;
    double dfm, pt_, s2pt, f, np;
    int i;

    dfm = lambda - x;
    pt_ = -log1pmx(dfm / x);
    s2pt = sqrt(2 * x * pt_);
    if (dfm < 0)
        s2pt = -s2pt;

    res12 = 0;
    res1_ig = res1_term = sqrt(x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++)
    {
        res12 += res1_ig * coefs_a[i];
        res12 += res2_ig * coefs_b[i];
        res1_term *= pt_ / i;
        res2_term *= 2 * pt_ / (2 * i + 1);
        res1_ig = res1_ig / x + res1_term;
        res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++)
    {
        elfb += elfb_term * coefs_b[i];
        elfb_term /= x;
    }
    if (!lower_tail)
        elfb = -elfb;
    f = res12 / elfb;

    if (lower_tail)
        np = 1. - pNorm(s2pt, 0.0, 1.0);
    else
        np = pNorm(s2pt, 0.0, 1.0);

    double nd = dNorm(s2pt, 0., 1.);

    return np + f * nd;
}

static double dpois_wrap(double x_plus_1, double lambda)
{
    static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;  // =3.196577e18

    if (!isfinite(lambda))
        return 0.;

    if (x_plus_1 > 1)
        return dpois_raw(x_plus_1 - 1, lambda);

    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
        return exp(-lambda - lgamma(x_plus_1));
    else
    {
        double d = dpois_raw(x_plus_1, lambda);
        return d * (x_plus_1 / lambda);
    }
}

double pgamma_raw(double x, double alph, int lower_tail)
{
    // Here, assume that  (x,alph) are not NA  &  alph > 0. 
    double res;

    if (x <= 0.)
        return 0.;
    if (x >= ML_POSINF)
        return 1.;

    if (x < 1)
    {
        res = pgamma_smallx(x, alph, lower_tail);
    }
    else if (x <= alph - 1 && x < 0.8 * (alph + 50))
    {
        double sum = pd_upper_series(x, alph);
        double d = dpois_wrap(alph, x);

        if (!lower_tail)
            res = 1 - d * sum;
        else
            res = sum * d;
    }
    else if (alph - 1 < x && alph < 0.8 * (x + 50))
    {
        double sum;
        double d = dpois_wrap(alph, x);

        if (alph < 1)
        {
            if (x * DBL_EPSILON > 1 - alph)
                sum = 1.;
            else
            {
                double f = pd_lower_cf(alph, x - (alph - 1)) * x / alph;
                sum = f;
            }
        }
        else
        {
            sum = pd_lower_series(x, alph - 1);
            sum = 1 + sum;
        }

        if (!lower_tail)
            res = sum * d;
        else
            res = 1 - d * sum;
    }
    else
    {
        res = ppois_asymp(alph - 1, x, !lower_tail);
    }

    if (res < DBL_MIN / DBL_EPSILON)
    {
        return exp(pgamma_raw(x, alph, lower_tail));
    }
    else
        return res;
}

double pgamma(double x, double alph, double scale, int lower_tail)
{

#ifdef IEEE_754
    if (isnan(x) || isnan(alph) || isnan(scale))
        return x + alph + scale;
#endif

    if (alph < 0. || scale <= 0.)
        return NAN;

    x /= scale;

#ifdef IEEE_754
    if (isnan(x)) // eg. original x = scale = +Inf
        return x;
#endif

    if (alph == 0.) // limit case.
        return (x <= 0) ? 0. : 1.;

    return pgamma_raw(x, alph, lower_tail);
}

#endif

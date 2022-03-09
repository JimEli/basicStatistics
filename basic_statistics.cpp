#define _USE_MATH_DEFINES 
#include <iostream>
#include <vector>

#include "common.h"
#include "normal.h"
#include "binomial.h"
#include "poisson.h"
#include "chisquare.h"
#include "student.h"

// Average for the data set.
template<typename T>
T mean(const std::vector<T>& v)
{
    if (v.empty())
        return { };

    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

//  Middle value of the data set.
template<typename T>
T median(std::vector<T> v)
{
    if (v.empty())
        return { };

    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin() + n, v.end());

    return v[n];
}

// Number that occurs the most in data set.
template<typename T>
T mode(std::vector<T> v)
{
    if (v.empty())
        return { };

    std::sort(v.begin(), v.end());

    // Finding max frequency.
    unsigned max_count = 1, count = 1;
    T r = *v.begin();

    for (auto it = v.begin(); it != v.end(); ++it)
    {
        if (it != v.begin() && *it == *(it - 1))
            count++;
        else
        {
            if (count > max_count)
            {
                max_count = count;
                r = *(it - 1);
            }
            count = 1;
        }
    }

    // When the last element is most frequent.
    if (count > max_count)
    {
        max_count = count;
        r = v[v.size() - 1];
    }

    // index = max_count;
    return r;
}

// Measure of how far a set of data are dispersed from their mean.
template<typename T>
T variance(const std::vector<T>& v)
{
    if (v.empty())
        return { };

    unsigned sz = v.size();
    const T mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
    auto variance_func = [&mean, &sz](T accumulator, const T& val) { return accumulator + ((val - mean) * (val - mean) / (sz - 1)); };

    return std::accumulate(v.begin(), v.end(), 0.0, variance_func);
}

// Measure of dispersement (tells how much data is spread out).
template<typename T>
T stdDeviation(std::vector<T> v)
{
    if (v.empty())
        return { };

    return sqrt(variance(v));
}

// Measure of how many standard deviations above/below the population mean.
template<typename T>
T zScore(const T x, const std::vector<T> v)
{
    if (v.empty())
        return { };

    return ((x - mean(v)) / stdDeviation(v));
}


// Sample usage.
void print(const std::string s, const double x) { std::cout << " " << s << " " << x << std::endl; }

int main()
{
    const std::vector<double> sample = { 3.0, 1.0, 5.0, 6.0, 3.0, 4.5 };

    std::cout << "mean: " << mean(sample);
    std::cout << ", median: " << median(sample);
    std::cout << ", mode:" << mode(sample) << std::endl;
    std::cout << "variance: " << variance(sample);
    std::cout << ", standard deviation: " << stdDeviation(sample) << std::endl;
    std::cout << "z-score:" << std::endl; std::for_each(sample.begin(), sample.end(), [sample](double x) { std::cout << "  " << x << ": " << zScore(x, sample) << std::endl; });


    // R Binomial Probabilities: P(X=k) = dbinom(#success, #trials, prob. success) and P(X<=k) = pbinom(#success, #trials, #prob. success)
    std::cout << "Binomial Distribution Probabilities\n";
    // A basketball player makes 44% of his/her 3-point shots. In a game he attempts 
    // ten 3-point shots. Find probability he makes exactly 7. = 0.067
    print("Exactly 7 shots:", dBinom(7, 10, 0.44));
    // 63% of adults drink coffee. Random sample of 25 adults is 
    // selected. Find probability more than 17 drink coffee. = 0.237
    // P(X>17) = P(X=18) + P(X=19) + ... P(X=25), as: dBinomial(18, 25, 0.63) + ... dBinomial(25, 25, 0.63)
    print("More than 17 drink coffee:", 1. - pBinom(17, 25, 0.63));


    // R Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Normal Distribution Probabilities\n";
    // Customers waiting time is normally distributed with mean of 2.58 minutes and standard deviation 
    // of 0.76. Find probability a random customer waits more than 4 minutes. = 0.969
    print("Waiting more than 4 mins:", pNorm(4, 2.58, 0.76));
    // Cholesterol levels are normally distributed with mean of 201 and standard 
    // deviation of 46. What cholesterol level separates the lowest 22%? = 165.488
    print("Normal distribution x:", qNorm(0.22, 201, 46));


    // Central Limit Theorem (means), n>=30 results in averages making approx. normal distribution with mu=mu, sigma=sigma/sqrt(n)
    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Central Limit Theorem (means)\n";
    // Mean cholesterol is 202 with standard deviation of 41. 110 person random 
    // sample what is probability for mean cholesterol between 190 and 200? = 0.303
    print("Cholesterol between 190-200:", pNorm(200, 202, qSigmaCLT(110, 41)) - pNorm(190, 202, qSigmaCLT(110, 41)));
    // Mean commute is 16 miles, with standard deviation of 8. Sample 75 commuters, 
    // there is 94% probability that mean commute is between what 2 distances? = 14.262 and 17.738
    std::cout << " 94% probability commute is between: " << qNorm((1. - 0.94) / 2, 16, qSigmaCLT(75, 8)); print(" to:", qNorm(1. - (1. - 0.94) / 2, 16, qSigmaCLT(75, 8)));


    // Central Limit Theorem (proportions), n*p>=5 and n*(1-p)>=5 can use normal distribution with mu=p, sigma=sqrt((p*(1-p))/n)
    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Central Limit Theorem (proportions)\n";
    // 63% of adults drink coffee daily. Random sample of 250 adults is selected. 
    // Find probability that more than 67% of sampled drink coffee daily. = 0.095
    print("Probability 67% of 250 drink coffee:", 1. - pNorm(0.67, 0.63, pSigmaCLT(250, 0.63)));
    // 68% of graduates have loan debt. From random sample of 85 grads, 
    // find probability that between 65 and 80% are in debt. = 0.713
    print("Probability 65-80% in debt:", pNorm(0.80, 0.68, pSigmaCLT(85, 0.68)) - pNorm(0.65, 0.68, pSigmaCLT(85, 0.68)));


    std::cout << "R Probability Normal Distribution (pnorm)\n";
    // In a precinct, 81% of the voters are registered Democrats. What is the probability that, in a random  
    // sample of 100 voters from this precinct, no more than 80 of the voters would be Democrats? = 0.3994
    print("No more than 80 out of 100 voters democrat:", pNormCDF(zPhat(100, 0.81, 0.8)));
    // Poll of 238 voters finds 137 approve an amendment. Construct a 95% confidence 
    // interval for proportion of voters supporting amendment. (0.513 to 0.638)
    double phat = (137. / 238.), MoE = E(238, Z95CI, phat); std::cout << "95% CI: " << phat - MoE << " to " << phat + MoE << " [margin of error: +/-" << MoE << "]" << std::endl;
    // "We are 95% confident that the 57.6% of the voters supporting the amendment is between 51.3% and 63.8%.


    std::cout << "Poisson Distributions\n";
    // Number of visits to a web page follows a Poisson distribution with mean 15 visits per hour.
    // What is probability of getting 10 or less visits per hour, P(X≤10)? =0.118
    // { double p=0.; for(int i=0; i<=10; i++) p+=dPois(i, 15); print("10 or less web hits per hour:", p); } // Equivalent
    print("10 or less web visits per hour:", pPois(10, 15));
    // What is probability of seeing exactly 3 blemishes on randomly selected piece of 
    // sheet metal, when on average one expects 1.2 blemishes, can be found with? =0.087
    print("3 blemishes on metal:", dPois(3, 1.2));


    std::cout << "dchisq=" << dchisq(5, 10);
    std::cout << ", pchisq=" << pchisq(5, 10);
    std::cout << ", qchisq=" << qchisq(.1, 200) << std::endl;


    std::cout << "qt=" << qt(.1, 5);
    std::cout << ", dt=" << dt(.1, 5);
    std::cout << ", pt=" << pt(.5, 100) << std::endl;
}

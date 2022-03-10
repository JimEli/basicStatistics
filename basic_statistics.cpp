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
    // 8.3% of all Americans have diabetes. Suppose that a random sample of 74 Americans is taken. What is the probability 
    // that the proportion of diabetics in the sample differs from the population proportion by less than 1%? = 0.2448
    print("Probability diabetics in sample of 74 differs from population by less than 1%:", 1. - 2 * pNormCDF(zPhat(74, 0.083, (0.083 - 0.01))));


    std::cout << "Confidence Intervals for Population Proportions\n";
    {
        {
            // Poll of 238 voters finds 137 approve an amendment. Construct a 95% confidence 
            // interval for proportion of voters supporting amendment. =(0.513 to 0.638)
            double phat = (137. / 238.), E = proportionMoE(238, Z95CI, phat);
            std::cout << " 95% CI: " << phat - E << " to " << phat + E << " [margin of error: +/-" << E << "]" << std::endl;
            // "We are 95% confident that the percentage of voters supporting the amendment is between 51.3% and 63.8% (57.6%).
        }
        {
            // Random sample of 250 college students determines their study hours have a sample mean of xbar=15.7 hours per week. 
            // If margin of error for the data using 95%CI is E=0.6 hours, construct a 95% confidence interval for the data. =(15.1 to 16.3)
            double phat = 15.7, E = 0.6;
            std::cout << " 95% CI: " << phat - E << " to " << phat + E << " [margin of error: +/-" << E << "]" << std::endl;
            // "We are 95% confident that the study hours are between 15.1 and 16.3 hrs (15.7)"
        }
        {
            // Survey to be conducteed of random sample asking whether they favor/oppose an amendment. How many should be polled 
            // to be sure that a 90% CI for the proportion who favor amendment will have a margin of error of 0.05.
            // A previous survey suggest the proportion in favor will be 84%. What sample size is needed? = 146
            double phat = 0.84, E = 0.05; std::cout << " 90% CI sample size: " << round(proportionN(E, Z90CI, phat) + .5) << std::endl;
            // Estimate sample size needed if no estimate of p is available. = 271
            phat = 0.5; std::cout << " 90% CI sample size (no prior estimate): " << round(proportionN(E, Z90CI, phat) + .5) << std::endl;
        }
        {
            // A newspaper wants to predict outcome of election by estimating proportion voters supporting particular 
            // candidate. What sample size is needed to yield an estimate within 3% with 97% CI? =1309
            double phat = 0.5, E = 0.03; std::cout << " 97% CI sample size: " << round(proportionN(E, qNorm(.97 + (1 - .97) / 2), phat) + .5) << std::endl;
        }
    }

    std::cout << "Confidence Intervals for Population Means\n";
    {
        {
            // Random survey sample of 43 renters found mean rent of $940 with 
            // standard deviation of $300. Construct a 91% CI for the mean rent. =(860.6 to 1019.4)
            unsigned n = 43 - 1; double xbar = 940, t = qt(1 - ((1 - .91) / 2), n), E = meanMoE(n, t, 300);
            std::cout << " 91% CI: " << xbar - E << " to " << xbar + E << " [margin of error: +/-" << E << "]" << std::endl;
            // "We are 91% confident that the mean rent is between $860.6 and $1019.4 ($940)"
        }
        {
            // We want an estimate of proportion of students who log into Facebook daily. 
            // Out of random sample of 200 students, 134 report logging in daily.
            // Construct 86% CI for proportion that log in daily. =(62.1 to 71.9)
            double phat = 134. / 200, E = proportionMoE(200, qNorm(.86 + (1 - .86) / 2), phat);
            std::cout << " 86% CI: " << phat - E << " to " << phat + E << " [margin of error: +/-" << E << "]" << std::endl;
            // "We are 86% confident that the mean log ins is between 62.1 and 71.9% (67%)"
        }
        {
            // to determine mean hours per week a person watches TV, how many people needed to estimate number hours within 2 hours 
            // with 94% CI? Results from previous survey indicate hours watched per week has a standard deviation of 7.5 hours. =50
            double z = qNorm(.94 + (1 - .94) / 2); std::cout << " 94% CI sample size: " << round(meanN(2, z, 7.5) + .5) << std::endl;
        }
    }


    std::cout << "Hypothesis Testing for Proportions\n";
    {
        {
            // McDonald's claims Monopoly game has 1 in 4 instant winners. 530 attempts, you win 112 (21%). Use alpha=0.05 significance.
            // State null and alternative hypothesis: H0 = "Game has win percentge of 25% (p=0.25)", H1 = "Game has win percentage < 25% (p<0.25)".
            // Find critical value and rejection region. =-1.645, left-tailed
            std::cout << " critcal z value: " << qNorm(0.05) << std::endl;
            // Find z statistic. = -2.073
            double phat = 112. / 530, p0 = 0.25; unsigned n = 530; double z = proportionHypothesisZ(n, phat, p0); print("z:", z);
            // Find p-value. =0.019
            double p = pNorm(z); print("p-value:", p);
            // Reject or accept H0? =[Our results (112 in 530) occur only 1.9% of time (.019<.05, or -2.073<-1.645), so reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.05);
            // State conclusion in sentence.
            std::cout << " At 0.05 level of significance, there is enough evidence to conclude that win rate is less than 25%.\n";
        }
    }


    std::cout << "Hypothesis Testing for Means\n";
    {
        {
            // Random sample of 50 pumpkins has mean circumference of 40.5cm with standard devication 1.6cm. Can you conclude 
            // mean circumference of all pumpkins is more than 40cm. Use alpha=0.01 significance level?
            // State null and alternative hypothesis: H0 = "Pumpkin mean circumference equals 40cm (mu=40)", H1 = "Pumpkin mean circumference is greater than 40cm (mu>40)".
            // Find critical value and rejection region. =2.405, right-tailed
            unsigned n = 50, DoF = n - 1; std::cout << " critcal t value: " << qt(1. - .01, DoF) << std::endl;
            // Find t statistic. =2.210
            double xbar = 40.5, mu = 40., sigma = 1.6, t = meanHypothesisT(n, xbar, mu, sigma); print("t:", t);
            // Find p-value. 0.016
            double p = 1. - pt(t, DoF); print("p-value:", p);
            // Reject or accept H0? =[Pumpkin mean circumferences >40 occur 1.6% of time (.016<.01, or 2.210<2.405), so do not reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.01);
            // State conclusion in sentence.
            std::cout << " At 0.01 level of significance, there is not enough evidence to conclude that all pumpkin mean circumference > 40cm.\n";

        }
    }
    
    
    std::cout << "Poisson Distributions\n";
    // Number of visits to a web page follows a Poisson distribution with mean 15 visits per hour.
    // What is probability of getting 10 or less visits per hour, P(X≤10)? =0.118
    // { double p=0.; for(int i=0; i<=10; i++) p+=dPois(i, 15); print("10 or less web hits per hour:", p); } // Equivalent
    print("10 or less web visits per hour:", pPois(10, 15));
    // What is probability of seeing exactly 3 blemishes on randomly selected piece of 
    // sheet metal, when on average one expects 1.2 blemishes, can be found with? =0.087
    print("3 blemishes on metal:", dPois(3, 1.2));


    std::cout << "Chi Square Distributions\n";
    std::cout << "dchisq=" << dchisq(5, 10);
    std::cout << ", pchisq=" << pchisq(5, 10);
    std::cout << ", qchisq=" << qchisq(.1, 200) << std::endl;
}

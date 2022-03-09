#define _USE_MATH_DEFINES 
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

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
    auto variance_func = [&mean, &sz](T accumulator, const T& val) { return accumulator + ((val - mean)*(val - mean)/(sz - 1)); };

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

    print("mean:", mean(sample));
    print("median:", median(sample));
    print("mode:", mode(sample));
    print("variance:", variance(sample));
    print("standard deviation:", stdDeviation(sample));
    std::cout << "z-score:" << std::endl; std::for_each(sample.begin(), sample.end(), [sample](double x) { std::cout << "  " << x << ": " << zScore(x, sample) << std::endl; });

    // R Binomial Probabilities: P(X=k) = dbinom(#success, #trials, prob. success) and P(X<=k) = pbinom(#success, #trials, #prob. success)
    std::cout << "Binomial Distribution Probabilities\n";
    // A basketball player makes 44% of his/her 3-point shots. In a game he attempts 
    // ten 3-point shots. Find probability he makes exactly 7. = 0.067
    print("Exactly 7 shots:", dBinom(7, 10, 0.44));
    // If he attempts ten 3-point shots, find probability he makes at least 7. P(x>=7) = P(x=7) + P(x=8) + P(x=9) + P(x=10).
    print("At least 7 shots:", dBinom(7, 10, 0.44) + dBinom(8, 10, 0.44) + dBinom(9, 10, 0.44) + dBinom(10, 10, 0.44));

    // 79% of adults need glasses. If 20 adults are randomly 
    // selected, find probability that 16 need glasses. = 0.217
    print("16 need glasses:", dBinom(16, 20, 0.79));

    // 63% of adults drink coffee. Random sample of 25 adults is 
    // selected. Find probability more than 17 drink coffee. = 0.237
    // P(X>17) = P(X=18) + P(X=19) + ... P(X=25), as: dBinomial(18, 25, 0.63) + ... dBinomial(25, 25, 0.63)
    print("More than 17 drink coffee:", 1. - pBinom(17, 25, 0.63));

    // Medical procedure produces side effects in 25% of patients. 60 people undergo the 
    // procedure. What is probability that fewer than 20 experience side effect? =0.908
    print("Fewer than 20 have side effects", pBinom(19, 60, 0.25));

    // An airline estimates 90% of pax show up for flight. If an airline sells 50 
    // tickets, what is the probability that not more than 46 show up? =0.750
    print("Not more than 46 show for airline seats:", pBinom(46, 50, 0.90));
    // An airline finds that 5% will not show up for flight. If the airline sells 160 tickets for a fight with 
    // 155 seats, what is the probability that a seat will be available for every person holding a reservation. =0.906
    print("All seats available:", pBinom(155, 160, 0.95));

    // Suppose 63% of people are right-handed. If 30 random people are selected, 
    // what is the probability that at least 15 of them is right-handed? =0.95
    print("At least 15 of 30 are right-handed:", 1. - pBinom(14, 30, 0.63));

    // 8% of vehicles fail emissions test. 12 random vehicles are 
    // tested. Find the probability that at least 1 fails test. =0.632
    print("At least 1 fails emissions:", 1. - dBinom(0, 12, 0.08));

    // McDonald's game produces 25% instant winners. Out of 530 
    // attempts, what is probability of only 108 winners? =0.007
    print("McDonald's winners:", pBinom(108, 530, 0.25));


    // R Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Normal Distribution Probabilities\n";
    // Customers waiting time is normally distributed with mean of 2.58 minutes and standard deviation 
    // of 0.76. Find probability a random customer waits more than 4 minutes. = 0.969
    print("Waiting more than 4 mins:", pNorm(4, 2.58, 0.76));

    // Cholesterol levels are normally distributed with mean of 201 and standard 
    // deviation of 46. What cholesterol level separates the lowest 22%? = 165.488
    print("Normal distribution x:", qNorm(0.22, 201, 46));

    // Survey of students study hours are normally distributed with a mean of 25 hours per week 
    // and standard deviation of 7. What proportion study more than 40 hours per week? = 0.016
    print("Study more than 40 hours:", 1. - pNorm(40, 25, 7));

    // Men's height normally distributed with mean of 68.6" and standard deviation 
    // of 2.8". Find height that separates the tallest 3% of men. = 73.867
    print("Tallest 3% of men:", qNorm((1. - 0.03), 68.6, 2.8));

    // A 12oz soda machine test shows it dispenses a normally distributed mean of 11.5oz with standard 
    // deviation of 0.21oz. What percent of cups will receive between 11.2 and 11.45ozs? = 0.329
    print("Soda machine dispenses 11.2 to 11.45oz:", pNorm(11.45, 11.5, 0.21) - pNorm(11.2, 11.5, 0.21));

    // Laptop batteries have mean life of 11.5 hours with standard deviation of 
    // 4.37 hours. Between what 2 lifetimes do 81% of batteries fall? = 5.771 and 17.229
    std::cout << " 81% of laptop batteries: " << qNorm((1. - .81) / 2, 11.5, 4.37) << " and " << qNorm(1. - (1. - .81) / 2, 11.5, 4.37) << std::endl;

    // Scores on a test were normally distributed with a mean of 75 and 
    // standard deviation of 8. Find 85th percentile of exams. = 83.288
    print("85th percentile of exams:", qNorm(0.85, 75, 8));


    // Central Limit Theorem (means), n>=30 results in averages making approx. normal distribution with mu=mu, sigma=sigma/sqrt(n)
    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Central Limit Theorem (means)\n";
    // Mean cholesterol is 202 with standard deviation of 41. 110 person random 
    // sample what is probability for mean cholesterol between 190 and 200? = 0.303
    print("Cholesterol between 190-200:", pNorm(200, 202, qSigmaCLT(110, 41)) - pNorm(190, 202, qSigmaCLT(110, 41)));

    // Typing speed is approx. normal, with mean speed of 45wpm and standard deviation of 10wpm. 
    // What is probability that random sample of 12 people have mean typing speed >40wpm? = 0.958
    print("12 people type greater than 40 wpm:", 1. - pNorm(40, 45, qSigmaCLT(12, 10)));

    // Mean commute is 16 miles, with standard deviation of 8. Sample 75 commuters, 
    // there is 94% probability that mean commute is between what 2 distances? = 14.262 and 17.738
    std::cout << " 94% probability commute is between: " << qNorm((1. - 0.94) / 2, 16, qSigmaCLT(75, 8)); print(" to:", qNorm(1. - (1. - 0.94) / 2, 16, qSigmaCLT(75, 8)));

    // Mean hourly wage is $14 with standard deviation of $8. A sample of 60 people is taken. 
    // There is a 32% probability that mean hourly wage in sample is greater than what? = 14.483
    print("32% probability mean hourly wage is greater than:", qNorm(1. - 0.32, 14, qSigmaCLT(60, 8)));

    // Credit scores are normally distributed with mean of 592 and a standard deviation of 
    // 106. What is the score of someone who is at the 92nd percentile? = 740.93
    print("92nd percentile credit score:", qNorm(0.92, 592, 106));


    // Central Limit Theorem (proportions), n*p>=5 and n*(1-p)>=5 can use normal distribution with mu=p, sigma=sqrt((p*(1-p))/n)
    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Central Limit Theorem (proportions)\n";
    // 63% of adults drink coffee daily. Random sample of 250 adults is selected. 
    // Find probability that more than 67% of sampled drink coffee daily. = 0.095
    print("Probability 67% of 250 drink coffee:", 1. - pNorm(0.67, 0.63, pSigmaCLT(250, 0.63)));

    // 68% of graduates have loan debt. From random sample of 85 grads, 
    // find probability that between 65 and 80% are in debt. = 0.713
    print("Probability 65-80% in debt:", pNorm(0.80, 0.68, pSigmaCLT(85, 0.68)) - pNorm(0.65, 0.68, pSigmaCLT(85, 0.68)));

    // Medical procedure produces side effects in 35% of patients. 80 people undergo the procedure. 
    // What is probability that fewer than 25 or fewer experience side effect? = 0.243
    print("Fewer than 25 have side effects:", pNorm((25. / 80.), 0.35, pSigmaCLT(80, 0.35)));


    std::cout << "Review\n";
    // SAT scores are normally distributed with mean of 514 and standard 
    // deviation of 118. What SAT score is required to be in top 2%? 
    print("SAT to be in top 2%:", qNorm(0.98, 514, qSigmaCLT(1, 118)));

    // 18% of commercials on TV are local advertisers. Of sample of 120 
    // commercials, what is probability that more than 20% are local? 
    print("Probability more than 20 local ads:", 1. - pNorm(0.20, 0.18, pSigmaCLT(120, 0.18)));

    // Lifetime of a tire is normally distriburted with mean of 40,000 miles and standard 
    // deviation of 5,000. Between what two lifetimes do 80% of lifetimes fall? 
    std::cout << " 80% of tire lifetime fall between: " << qNorm(0.10, 40000, qSigmaCLT(1, 5000)); print(" to:", qNorm(0.90, 40000, qSigmaCLT(1, 5000)));

    // Lightbulb has mean lifetime of 1500 hours with standard deviation = 100 hours. Out of 
    // sample of 50, what is probability mean lifetime is between 1490 and 1550 hours? 
    print("Probability bulb lifetime is between 1490-1550hrs:", pNorm(1550, 1500, qSigmaCLT(50, 1500)) - pNorm(1490, 1500, qSigmaCLT(50, 1500)));

    // An airline estimates 90% of pax show up for a flight. If airline sells 
    // 105 tickets, what is the probability 100 or fewer pax show up? 
    print("Fewer than 100 show for Airline seats:", pNorm((100. / 105.), 0.90, pSigmaCLT(105, 0.90)));


    std::cout << "R Probability Normal Distribution (pnorm)\n";
    // In a precinct, 81% of the voters are registered Democrats. What is the probability that, in a random  
    // sample of 100 voters from this precinct, no more than 80 of the voters would be Democrats? = 0.3994
    print("No more than 80 out of 100 voters democrat:", pNormCDF(zPhat(100, 0.81, 0.8)));
    // In a different precinct, 79% of the voters are registered Republicans. What is the probability that in a 
    // random sample of 100 voters from this precinct at least 68 of the voters would be Republicans? = 0.9965 
    print("At least 68 out of 100 veters republican:", 1. - pNormCDF(zPhat(100, 0.79, 0.68)));
    // 8.3% of all Americans have diabetes. Suppose that a random sample of 74 Americans is taken. What is the probability 
    // that the proportion of diabetics in the sample differs from the population proportion by less than 1%? = 0.2448
    print("Probability diabetics in sample of 74 differs from population by less than 1%:", 1. - 2 * pNormCDF(zPhat(74, 0.083, (0.083 - 0.01))));
    // 8.3% of all Americans have diabetes. Suppose that a random sample of 91 Americans is taken. What is the probability that 
    // the proportion of people in the sample who are diabetic differs from the population proportion by more than 2%? = 0.4892
    print("Probability diabetics in sample of 91 differs from population by more than 2%:", 2 * pNormCDF(zPhat(91, 0.083, (0.083 - 0.02))));



    {
        // Mean cholesterol is 197mg/dL, with a standard deviation of 35 and levels are normally distributed. 
        // What is probability that a randomly selected sample of 150 adults will have a mean of less than 183?
        //double z = (183 - 197) / (35 / sqrt(150)), p = pNorm(183, 197, 35); std::cout << z << ", " << p << ", " << std::endl;
    }
    {
        // Poll of 238 voters finds 137 approve an amendment. Construct a 95% confidence 
        // interval for proportion of voters supporting amendment. (0.513 to 0.638)
        double phat = (137. / 238.), MoE = E(238, Z95CI, phat); std::cout << "95% CI: " << phat - MoE << " to " << phat + MoE << " [margin of error: +/-" << MoE << "]" << std::endl;
        // "We are 95% confident that the 57.6% of the voters supporting the amendment is between 51.3% and 63.8%.
    }
    // A random sample of 250 college students determines that their study hours have a sample mean of xbar=15.7 hours per week. 
    // If the margin of error for her data using a 95%CI is E=0.6 hours, construct a 95% confidence interval for her data. ()
    //phat = 15.7, MoE = E(250, Z95CI, phat); //std::cout << "95% CI: " << phat - MoE << " to " << phat + MoE << " [margin of error: +/-" << MoE << "]" << std::endl;
    //double E(n, z, phat) { return (z * sqrt(phat * (1. - phat) / n)); } 




    // Survey to be conducteed of random sample asking whether they favor/oppose an amendment. How many should be polled 
    // to be sure that a 90% CI for the proportion who favor amendment will have a margin of error of 0.05? =

    // A previous survey suggest the proportion in favor will be 84 % .What smaple size is needed? =

    // Estimate sample size needed if no estimate of p is available. =




    std::cout << "Poisson Distributions\n";
    // Number of visits to a web page follows a Poisson distribution with mean 15 visits per hour.
    // What is probability of getting 10 or less visits per hour, P(X≤10)? =0.118
    // { double p=0.; for(int i=0; i<=10; i++) p+=dPois(i, 15); print("10 or less web hits per hour:", p); } // Equivalent
    print("10 or less web visits per hour:", pPois(10, 15));

    // A call center gets average 13 calls per hour. What is probability that 
    // in 15 minutes call center will receive less than 6 call? =0.889
    print("Call center gets less than 6 calls:", pPois(5, 13. / 4));

    // What is probability of seeing exactly 3 blemishes on randomly selected piece of 
    // sheet metal, when on average one expects 1.2 blemishes, can be found with? =0.087
    print("3 blemishes on metal:", dPois(3, 1.2));

    // 170 oil tankers arrived randomly and independently at port over last 104 days. 
    // Only 2 tankers may be unloaded per day, extra oil tankers wait.
    // What is the probability that no tankers will arrive on Tuesday? =0.195
    // What is the probability that more than two will arrive on Friday, 
    // so that some will wait for Saturday to be unloaded? =0.225
    // Assuming no tankers are left over from Tuesday, what is the probability that exactly one 
    // tanker will be left over from Wednesday and none will be left over from Thursday? =0.142
    print("No tankers on Tuesday:", pPois(0, 170. / 104));
    print("2 tankers on Friday:", 1. - pPois(2, 170. / 104));
    print("3 tankers on Wednesday:", dPois(3, 170. / 104));


    /*
    union
    {
        unsigned int ui[2];
        double d;
    } ud;

    ud.ui[0] = 0x7fffffff;
    ud.ui[1] = 0xffffffff;

    std::cout << "float value is: " << ud.d << std::endl;
    std::cout << std::hex << "float value is: " << ud.ui[0] << ":" << ud.ui[1] << std::dec << std::endl;
  */



  //  std::cout << "\n" << dNorm(10, 10, 1) << std::endl;  // density f(x) [0.3989]
  //  std::cout << pNorm(10, 10, 1) << std::endl;  // probability at x=10 [50%]
  //  std::cout << qNorm(.5, 10, 1) << std::endl;  // x at probaility=0.5 [10]
  //std::cout << dBinom(10, 20, 0.15) << std::endl;
  //std::cout << pBinom(10, 20, 0.85) << std::endl;
  //std::cout << dPois(2, 10) << std::endl; // 0.002269
  //std::cout << pPois(2, 10) << std::endl; // 0.002769

  //  std::cout << dchisq(5, 10) << std::endl;
  //  std::cout << pchisq(5, 10) << std::endl;
  //  std::cout << qchisq(.1, 200) << std::endl;

  std::cout << "qt=" << qt(.1, 5) << std::endl;
  std::cout << "dt=" << dt(.1, 5) << std::endl;
  std::cout << "pt=" << pt(.5, 100) << std::endl;

}

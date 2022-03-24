#define _USE_MATH_DEFINES 
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
// Our stat package includes.
#include "common.h"
#include "normal.h"
#include "binomial.h"
#include "student.h"
#include "chisquare.h"
#include "poisson.h"

// Sample usage.
void print(const std::string& s, const double x) { std::cout << " " << s << " " << x << std::endl; }
void printCI(const std::string ci, const double phat, const double E) { std::cout << " " << ci << "% CI: " << phat - E << " to " << phat + E << " [margin of error: +/-" << E << "]\n"; }

int main()
{
    std::vector<double> sample = { 3.0, 1.0, 5.0, 6.0, 3.0, 4.5 };

    print("mean:", mean(sample));
    print("median:", median(sample));
    print("mode:", mode(sample));
    print("variance:", variance(sample));
    print("standard deviation:", standardDeviation(sample));
    std::cout << "z-score:" << std::endl; std::for_each(sample.begin(), sample.end(), [&sample](double x) { std::cout << "  " << x << ": " << zScore(x, sample) << std::endl; });


    // Binomial Probabilities: P(X=k) = dbinom(#success, #trials, prob. success) and P(X<=k) = pbinom(#success, #trials, #prob. success)
    std::cout << "Binomial Distribution Probabilities\n";
    {
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
    }

    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Normal Distribution Probabilities\n";
    {
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
    }

    // Central Limit Theorem (means), n>=30 results in averages making approx. normal distribution with mu=mu, sigma=sigma/sqrt(n)
    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Central Limit Theorem (means)\n";
    {
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
    }

    // Central Limit Theorem (proportions), n*p>=5 and n*(1-p)>=5 can use normal distribution with mu=p, sigma=sqrt((p*(1-p))/n)
    // Normal Distribution Probabilities (percent, proportion): pnorm(z) x->z->area, qnorm(left area) area->z->x
    std::cout << "Central Limit Theorem (proportions)\n";
    {
        // 63% of adults drink coffee daily. Random sample of 250 adults is selected. 
        // Find probability that more than 67% of sampled drink coffee daily. = 0.095
        print("Probability 67% of 250 drink coffee:", 1. - pNorm(0.67, 0.63, pSigmaCLT(250, 0.63)));

        // 68% of graduates have loan debt. From random sample of 85 grads, 
        // find probability that between 65 and 80% are in debt. = 0.713
        print("Probability 65-80% in debt:", pNorm(0.80, 0.68, pSigmaCLT(85, 0.68)) - pNorm(0.65, 0.68, pSigmaCLT(85, 0.68)));

        // Medical procedure produces side effects in 35% of patients. 80 people undergo the procedure. 
        // What is probability that fewer than 25 or fewer experience side effect? = 0.243
        print("Fewer than 25 have side effects:", pNorm((25. / 80.), 0.35, pSigmaCLT(80, 0.35)));
    }

    std::cout << "Review\n";
    {
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
    }

    std::cout << "R Probability Normal Distribution (pnorm)\n";
    {
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
    }

    std::cout << "Confidence Intervals for 1-Sample Population Proportions\n";
    {
        {
            // Poll of 238 voters finds 137 approve an amendment. Construct a 95% confidence 
            // interval for proportion of voters supporting amendment. =(0.513 to 0.638)
            printCI("95", (137. / 238.), proportionMoE(238, Z95CI, (137./238.)));
            // "We are 95% confident that the percentage of voters supporting the amendment is between 51.3% and 63.8% (57.6%).
        }
        {
            // Random sample of 250 college students determines their study hours have a sample mean of xbar=15.7 hours per week. 
            // If margin of error for the data using 95%CI is E=0.6 hours, construct a 95% confidence interval for the data. =(15.1 to 16.3)
            double phat = 15.7, E = 0.6; printCI("95", phat, E);
            // "We are 95% confident that the study hours are between 15.1 and 16.3 hrs (15.7)"
        }
        {
            // Survey to be conducteed of random sample asking whether they favor/oppose an amendment. How many should be polled 
            // to be sure that a 90% CI for the proportion who favor amendment will have a margin of error of 0.05.
            // A previous survey suggest the proportion in favor will be 84%. What sample size is needed? = 146
            double phat = 0.84, E = 0.05; printCI("90", phat, round(proportionN(E, Z90CI, phat) + .5));
            // Estimate sample size needed if no estimate of p is available. = 271
            phat = 0.5; printCI("90", phat, round(proportionN(E, Z90CI, phat) + .5));
        }
        {
            // A newspaper wants to predict outcome of election by estimating proportion voters supporting particular 
            // candidate. What sample size is needed to yield an estimate within 3% with 97% CI? =1309
            double phat = 0.5, E = 0.03; printCI("97", phat, round(proportionN(E, qNorm(.97 + (1 - .97) / 2), phat) + .5));
        }
    }

    std::cout << "Confidence Intervals for 1-Sample Population Means\n";
    {
        {
            // Random survey sample of 43 renters found mean rent of $940 with 
            // standard deviation of $300. Construct a 91% CI for the mean rent. =(860.6 to 1019.4)
            unsigned n = 43 - 1; double xbar = 940, t = qt(1 - ((1 - .91) / 2), n), E = meanMoE(n, t, 300); printCI("91", xbar, E);
            // "We are 91% confident that the mean rent is between $860.6 and $1019.4 ($940)"
        }
        {
            // We want an estimate of proportion of students who log into Facebook daily. 
            // Out of random sample of 200 students, 134 report logging in daily.
            // Construct 86% CI for proportion that log in daily. =(62.1 to 71.9)
            double phat = 134. / 200, E = proportionMoE(200, qNorm(.86 + (1 - .86) / 2), phat); printCI("86", phat, E);
            // "We are 86% confident that the mean log ins is between 62.1 and 71.9% (67%)"
        }
        {
            // to determine mean hours per week a person watches TV, how many people needed to estimate number hours within 2 hours 
            // with 94% CI? Results from previous survey indicate hours watched per week has a standard deviation of 7.5 hours. =50
            double z = qNorm(.94 + (1 - .94) / 2); printCI("94", z, round(meanN(2, z, 7.5) + .5));
        }
    }


    // Type I Error: When we reject H0, but in reality H0 is true. 
    // * Probability of this error is alpha. 
    // * Reduce probability of this error by lowering alpha.
    // Type II Error: When we don't reject H0, but in reality H0 was false. 
    // * Probability of this error is beta. 
    // * Reduce probability of this error by increasing sample size.
    std::cout << "1-Sample Hypothesis Testing for Proportions\n";
    {
        {
            // McDonald's claims Monopoly game has 1 in 4 instant winners. 530 attempts, you win 112 (21%). Use alpha=0.05 significance.
            // State null and alternative hypothesis: H0 = "Game has win percentge of 25% (p=0.25)", H1 = "Game has win percentage < 25% (p<0.25)".
            // Find critical value and rejection region. =-1.645, left-tailed
            std::cout << " critcal z value: " << qNorm(0.05) << std::endl;
            // Find z statistic. = -2.073
            unsigned n = 530; double phat = 112. / n, p0 = 0.25, z = proportionHypothesisZ(n, phat, p0); print("z:", z);
            // Find p-value. =0.019
            double p = pNorm(z); print("p-value:", p);
            // Reject or accept H0? =[Our results (112 in 530) occur only 1.9% of time (.019<.05, or -2.073<-1.645), so reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.05);
            // State conclusion in sentence.
            std::cout << " At 0.05 level of significance, there is enough evidence to conclude that win rate is less than 25%.\n";
        }
        {
            // Survey reports 15% of college student read newspaper daily. A researcher wants to know if percentage of newspaper readers among 
            // students at Specific college differs from percentage in general. A random sample survey of 200 at Specific college finds 40 read 
            // newspaper daily. Can it be concluded that proportion of reders at Specific college differs from 15%? Use alpha=0.03 level of significance.
            // State null and alternative hypothesis: H0 = "Specific college newspaper readers equal 15% (p=0.15)", H1 = "Specfic college newspaper readers does not equal 15% (p!=0.15)".
            // Find critical value and rejection region. =+/-2.170, two-tailed (right & left)
            std::cout << " critcal z value: " << qNorm(.03 / 2) << " and " << qNorm(1. - (.03 / 2)) << std::endl;
            // Find z statistic. =-1.980
            unsigned n = 200; double phat = 40. / n, p0 = 0.15, z = proportionHypothesisZ(n, phat, p0); print("z:", z);
            // Find p-value. =0.048
            double p = 2. * pNorm(-z); print("p-value:", p);
            // Reject or accept H0? =[Results (40 in 200) occur only 4.8% of time (0.048<0.03, or -1.98<-2.17), so don't reject H0]
            std::cout << " We ";  DecideHypothesis(p, .03);
            // State conclusion in sentence.
            std::cout << " At 0.03 level of significance, there is not enough evidence to conclude that Specific college students daily newspaper readers differs from 15%.\n";
        }
    }


    std::cout << "1-Sample Hypothesis Testing for Means\n";
    {
        {
            // Random sample of 50 pumpkins has mean circumference of 40.5cm with standard devication 1.6cm. Can you conclude 
            // mean circumference of all pumpkins is more than 40cm. Use alpha=0.01 significance level?
            // State null and alternative hypothesis: H0 = "Pumpkin mean circumference equals 40cm (mu=40)", H1 = "Pumpkin mean circumference is greater than 40cm (mu>40)".
            // Find critical value and rejection region. =2.405, right-tailed
            unsigned n = 50, DoF = n - 1; std::cout << " critcal t value: " << qt(1. - .01, DoF) << std::endl;
            // Find t statistic. =2.210
            double xbar = 40.5, mu = 40., sigma = 1.6, t = meanHypothesisT(n, xbar, mu, sigma); print("t:", t);
            // Find p-value. =0.016
            double p = 1. - pt(t, DoF); print("p-value:", p);
            // Reject or accept H0? =[Pumpkin mean circumferences >40 occur 1.6% of time (.016<.01, or 2.210<2.405), so do not reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.01);
            // State conclusion in sentence.
            std::cout << " At 0.01 level of significance, there is not enough evidence to conclude that all pumpkin mean circumference > 40cm.\n";
        }
        {
            // Mean annual tuition for a sample of 12 colleges was $39,200 with standard deviation of $4300. 
            // Can you conclude that mean tuition for private colleges is differnet from $35,500. Use alpha=0.03 significance level.
            // State null and alternative hypothesis: H0 = "Private college mean tuition equals $35,500 (mu=35500)", H1 = "Private college mean tuition does not equal $35,500 (mu!=35500)".
            // Find critical value and rejection region. = +/-2.491, two-tailed (left & right)
            unsigned n = 12, DoF = n - 1; std::cout << " critcal t value: " << qt(.03 / 2, DoF) << " and " << qt(1. - (.03 / 2), DoF) << std::endl;
            // Find t statistic. =-2.981
            double xbar = 39200., mu = 35500., sigma = 4300., t = meanHypothesisT(n, xbar, mu, sigma); print("t:", t);
            // Find p-value. =0.013
            double p = 2. * pt(-t, DoF); print("p-value:", p);
            // Reject or accept H0? =[Mean college tuituion equals 35500 at 0.6% of time (.013<.03, or -2.98<-2.491), so do not reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.03);
            // State conclusion in sentence.
            std::cout << " At 0.03 level of significance, is enough evidence to conclude that college mean tuituion is different from $35,500.\n";
        }
    }


    std::cout << "Two-Sample Hypothesis Testing for Proportions\n";
    {
        {
            // A new medicine test with sample of 2529 people, 1304 in control group given placebo and 1225 in treatment group 
            // given medicine. 150 in control group and 113 in treatment group got sick. Can company conclude that proportion 
            // who got sick differs between control/treatment groups? Use alpha=0.07 level of significance. Control group=p1.
            // State null and alternative hypothesis: H0 = "Proportion who got sick in placebo and treatment groups are equal (p1=p2)", H1 = "... placebo and teratment groups are not equal (p1!=p2)".
            // Find critical value and rejection region. =+/-1.812, two-tailed (left & right)
            std::cout << " critcal z value: " << qNorm(0.07 / 2) << " and " << qNorm(1. - (0.07 / 2)) << std::endl;
            // Find z statistic. = 1.894
            unsigned n1 = 1304, n2 = 1225; double phat1 = 150. / n1, phat2 = 113. / n2, z = proportionHypothesisZ2(n1, n2, phat1, phat2); print("z:", z);
            // Find p-value. = 0.058
            double p = 2. * pNorm(-z); print("p-value:", p);
            // Reject or accept H0? =[Differnece between group could have occurred 5.8% of time (.061<.07, or -1.894<-1.812), so reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.07);
            // State conclusion in sentence.
            std::cout << " At 0.07 level of significance, there is enough evidence to conclude that proportion who got sick in Placebo differs from treatment group.\n";
            // What type of error could have occurred, and what is the probablility? =Type I Error, probability=0.07 (alpha)
        }
    }

    std::cout << "Two-Sample Hypothesis Testing for Means\n";
    {
        {
            // Random sample of 17 business majors had mean 2.81 GPA with standard deviation of 0.27. A random survey of 23 
            // psychology majors had mean 2.97 GPA with standard deviation of 0.23. Can you conclude that mean GPA of psychology 
            // majors is greater than business majors? Use alpha = 0.02 significance level. Psychology majors=mu1. 
            // State null and alternative hypothesis: H0 = "Mean GPA of business and psychology majors are equal (mu1=mu2)", H1 = "Mean GPA of psychology majors is greater than business majors (mu1>mu2)".
            // Find critical value and rejection region. =2.235, right-tailed
            unsigned n1 = 23, n2 = 17, DoF = n2 - 1; std::cout << " critcal t value: " << qt(1. - .02, DoF) << std::endl;
            // Find t statistic. =1.971
            double xbar1 = 2.97, xbar2 = 2.81, sigma1 = 0.23, sigma2 = 0.27, t = meanHypothesisT2(n1, n2, xbar1, xbar2, sigma1, sigma2); print("t:", t);
            // Find p-value. =0.033
            double p = 1. - pt(t, DoF); print("p-value:", p);
            // Reject or accept H0? =[Mean psychology GPA > business 3.3% of time (.033<.02, or 1.971>2.235), so don't reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.02);
            // State conclusion in sentence.
            std::cout << " At 0.02 level of significance, there is not enough evidence to conclude mean psychology GPAs > business GPAs.\n";
            // What type of error could have occurred, and what is the probablility? =Type II Error, probability=xx (beta)
        }
        {
            // Concetration of benzene in 5 random samples of untreated wastewater had mean of 7.8mg/L with standard deviation 1.4mg/L. 
            // 7 random specimens of treated wastewater had means of 3.2mg/L with standard deviation 1.7mg/L. Assume both samples from 
            // approx normal populations. Can you conclude mean concentration is less in treated wastewater vs untreated? Use alpha=0.03 
            // level. mu1=treated, m2=untreated. 
            // State null and alternative hypothesis: H0 = "Mean benzene in untreate equals treated (mu1=mu2)", H1 = "Mean concentration of treated is less than untreated (mu1>mu2)".
            // Find critical value and rejection region. =-2.601, left-tailed
            unsigned n1 = 7, n2 = 5, DoF = n2 - 1; std::cout << " critcal t value: " << qt(.03, DoF) << std::endl;
            // Find t statistic. =-5.127
            double xbar1 = 3.2, xbar2 = 7.8, sigma1 = 1.7, sigma2 = 1.4, t = meanHypothesisT2(n1, n2, xbar1, xbar2, sigma1, sigma2); print("t:", t);
            // Find p-value. =0.003
            double p = pt(t, DoF); print("p-value:", p);
            // Reject or accept H0? =[Mean psychology GPA > business 3.3% of time (.033<.02, or 1.971>2.235), so don't reject H0]
            std::cout << " We ";  DecideHypothesis(p, 0.02);
            // State conclusion in sentence.
            std::cout << " At 0.03 level of significance, there is enough evidence to conclude treated wastewater benzene is less than untreated.\n";
            // What type of error could have occurred, and what is the probablility? =Type I Error, probability=alpha (0.03)
        }
    }

    std::cout << "Two-Sample (Independent) Hypothesis Testing for Means\n";
    {
        {
            // Do business travelers walk at different speed than leisure travelers? Researcher measured 12 of each (fpm). 
            std::vector<double> bt = { 272, 267, 249, 240, 382, 303, 289, 253, 327, 273, 216, 326 };
            std::vector<double> lt = { 328, 249, 255, 242, 210, 243, 212, 263, 254, 265, 270, 261 };
            // Can the researcher conclude business travelers walk different speed than leisure travelers? Use alpha=0.9 level. mu1=business
            // State null and alternative hypothesis: H0 = "Business and leisure travelers walk same speed (mu1=mu2)", H1 = "Business and leisure travelers walk at different speed (mu1!=mu2)".
            // Find critical value and rejection region. =+/-1.859, two-tailed
            unsigned n1 = bt.size(), n2 = lt.size(), DoF = /*smaller sample*/n2 - 1; std::cout << " critcal t value: " << qt(.09 / 2, DoF) << " and " << qt(1. - (.09 / 2), DoF) << std::endl;
            // Find t statistic. =1.822
            double xbar1 = mean(bt), xbar2 = mean(lt), sigma1 = standardDeviation(bt), sigma2 = standardDeviation(lt);
            double t = meanHypothesisT2(n1, n2, xbar1, xbar2, sigma1, sigma2);
            print("t:", t);
            // Find p-value. =0.096
            double p = 2. * pt(-t, DoF);
            print("p-value:", p);
            // Reject or accept H0? =don't reject
            std::cout << " We ";  DecideHypothesis(p, 0.02);
            // State conclusion in sentence:
            std::cout << " At 0.09 level of significance, there is not enough evidence to conclude mean speeds differ between travelers.\n";
        }
    }


    std::cout << "Two-Sample (Dependent/Matched-Pair) Hypothesis Testing for One Mean\n";
    {
        {
            // Do classrooms contain more bacteria in spring or fall? 8 classrooms randomly selected and measure amount of bacteria in smae rooms spring and fall.
            std::vector<double> fall = { 8.4, 11.3, 6.9, 14.8, 14.5, 11.7, 12.8, 9.7 };
            std::vector<double> spring = { 8.1, 12.5, 3.4, 10.2, 13.3, 12.2, 7.3, 10.8 }, delta;
            // Can we conclude classrooms have less bacteria in spring than fall? Use alpha = 0.08 level of significance. mu=differences: fall>spring=0, fall-spring>0
            for (size_t i = 0; i < fall.size(); i++)
                delta.push_back(fall[i] - spring[i]);
            unsigned n = fall.size(); double mu = 0., xbar = mean(delta), sigma = standardDeviation(delta);
            // State null and alternative hypothesis: H0 = "Spring and fall classrooms have same amount of bacteria (mu=0)", H1 = "Classroom have less bacteria in spring (mu>0)".
            // Find critical value and rejection region. =1.571
            unsigned DoF = n - 1; std::cout << " critcal t value: " << qt(1. - .08, DoF) << std::endl;
            // Find t statistic. =1.641
            double t = meanHypothesisT(n, xbar, mu, sigma); print("t:", t);
            // Find p-value. 
            double p = 1. - pt(t, DoF); print("p-value:", p);
            // Reject or accept H0? =reject H0
            std::cout << " We ";  DecideHypothesis(p, 0.08);
            // State conclusion in sentence.
            std::cout << " At 0.08 level of significance, there is enough evidence to conclude spring classrooms have less bacteria than fall.\n";
        }
    }


    std::cout << "Confidence Intervals for 2-Sample Population Proportions\n";
    {
        {
            // Study conducted comparing effectieness of distance vs. classroom learning. 12 students took business admin course online, 14 took in classroom. Final exams are below.
            std::vector<double> online = { 66, 75, 85, 64, 88, 77, 74, 91, 72, 69, 77, 83 };
            std::vector<double> classroom = { 80, 83, 64, 81, 75, 80, 86, 81, 91, 64, 99, 85, 74, 77 };
            // Construct 95% confidence interval for difference bewteen mean scores.
            unsigned n1 = online.size(), n2 = classroom.size(); double xbar1 = mean(online), xbar2 = mean(classroom), sigma1 = standardDeviation(online), sigma2 = standardDeviation(classroom);
            // 1. Find critical area. =+/-2.201
            unsigned DoF = std::min(n1, n2) - 1; double t = qt(1. - (1 - .95) / 2, DoF); std::cout << " t: " << qt((1 - .95) / 2, DoF) << " and " << t << std::endl;
            // 2. Margin of error =7.740
            double E = meanMoE2(n1, n2, sigma1, sigma2, t); print("margin of error:", E);
            // 3. MoE
            printCI("95", xbar1 - xbar2, E);
            // " We are 95% confident that the mean scores for the 2 types of instruction is between -10.99 and 4.49"
        }
        {
            // A simple random sample of 90 store receipts a week before an ad campaign found 21 pretzel purchases. 
            // Another sample a week after the ad campaign of 70 receipts had 39 purchases. 
            // Construct a 98% CI for the difference between proportions of customers purchasing pretzels before and after ad campaign. 
            unsigned n1 = 90, n2 = 70; double phat1 = 21. / 90, phat2 = 39. / 70, phat = phat2 - phat1 /* after - before */;
            // 1. critical area =+/-2.326
            double z = qNorm((1. - 0.98) / 2); std::cout << " z: " << z << " and " << -z << std::endl;
            // 2. Margin of error =0.173
            double E = proportionMoE2(n1, n2, -z, phat1, phat2); print("margin of error:", E); printCI("98", phat, E);
            // "We are 98% confident that the difference between proportion of pretzels purchased before and after ad campagain is between 15.1% and 49.7%"
        }
        {
            // 8 individuals with high cholesterol took a new drug. Before and after measurements below:
            std::vector<double> before = { 283, 299, 274, 284, 248, 275, 293, 277 };
            std::vector<double> after = { 215, 206, 187, 212, 178, 212, 192, 196 }, delta;
            // Construct a 90% CI for the mean reduction.
            for (size_t i = 0; i < before.size(); i++)
                delta.push_back(before[i] - after[i]);
            unsigned n = after.size(), DoF = n - 1; double xbar = mean(delta), sigma = standardDeviation(delta);
            // 1. Find critical values. =+/-1.895
            double t = qt((1 - .90) / 2, DoF); std::cout << " t: " << t << " and " << qt(1 - ((1 - .90) / 2), n) << std::endl;
            // 2. Find margin of error. =8.967
            double E = meanMoE(n, -t, sigma); print("margin of error:", E); printCI("90%", xbar, E);
            // "We are 90% confident that the mean reduction in cholesterol is between 70.408 and 88.342."
        }
        {
            // Survey of college students reported out of 413 male the average energy drinks consumed per month was 2.49 
            // with a standard deviation of 4.87. Out of 382 female average was 1.22 with standard deviation of 3.24. 
            // Construct 99% CI for difference between men and women in mean energy drinks consumed.
            unsigned n1 = 413, n2 = 382; double xbar1 = 2.49, xbar2 = 1.22, sigma1 = 4.87, sigma2 = 3.24;
            // 1. Find critical area. =+/-
            unsigned DoF = std::min(n1, n2) - 1; double t = qt(1. - (1 - .95) / 2, DoF); std::cout << " t: " << qt((1 - .95) / 2, DoF) << " and " << t << std::endl;
            // 2. Margin of error =
            double E = meanMoE2(n1, n2, sigma1, sigma2, t); print("margin of error:", E);
            // 3. MoE
            printCI("99", xbar1 - xbar2, E);
            // " We are 99% confident that the difference between men and women mean energy drinks consumed is between 0.697 and 1.843"
        }
    }


    std::cout << "Hypothesis Test of Goodness of Fit (chi-square)\n";
    {
        {
            // Biology professor claims on average 15% of students get A, 30% B, 40% C, 10% D and 5% F. 
            // Grades of random sample in following order: A, B, C, D, F
            std::vector<double> grades = { .15, .30, .40, .10, .05 }, observed = { 8., 35., 49., 15., 11. }, expected;
            // Can you conclude grades follow different distribution from claimed? Use alpha = 0.02 significance level.
            // State null and alternative hypothesis: H0 = "observed = expected", H1 = "observed != expected".
            // Find critical value, always rigth-tailed! =11.668
            unsigned DoF = grades.size() - 1; double cv = qchisq(1. - 0.02, DoF); print("critical value: ", cv);
            // Find chi-square. =10.665
            double n = std::accumulate(observed.begin(), observed.end(), 0.);
            for (size_t i = 0; i < grades.size(); i++)
                expected.push_back(grades[i] * n);
            double cs = findChiSquare<double>(observed, expected); print("chi-square: ", cs);
            // Find p-value. =0.031
            double p = 1. - pchisq(cs, DoF); print("p-value: ", p);
            // Reject or accept H0? =reject H0
            std::cout << " We ";  DecideHypothesis(p, 0.02);
            // State conclusion in sentence.
            std::cout << " At 0.02 level of significance, there is not enough evidence to conclude grades follow distribution different than claimed.\n";
        }
        {
            // Total absences on each day of week (Mon thru Fri) for a class.
            std::vector<double> observed = { 39, 24, 20, 26, 41 }, expected(observed.size()), csq;
            // Is there evidence to conclude absences are more likely on some days than others? Use alpha=0.02 significance level.
            // State null and alternative hypothesis: H0 = "observed = expected", H1 = "observed != expected".
            // Find critical value, always rigth-tailed! =11.668
            unsigned DoF = observed.size() - 1; double cv = qchisq(1. - 0.02, DoF); print("critical value: ", cv);
            // Find chi-square. =11.8
            double n = std::accumulate(observed.begin(), observed.end(), 0.);
            std::fill(expected.begin(), expected.end(), n / 5);
            double cs = findChiSquare<double>(observed, expected); print("chi-square: ", cs);
            // Find p-value. =0.019
            double p = 1. - pchisq(cs, DoF); print("p-value: ", p);
            // Reject or accept H0? =reject H0
            std::cout << " We ";  DecideHypothesis(p, 0.02);
            // State conclusion in sentence.
            std::cout << " At 0.02 level of significance, there is enough evidence to conclude absences are more likely on some days than others.\n";
        }
    }


    std::cout << "Hypothesis Test for Independence (chi-square)\n";
    {
        {
            // Table below is a study to determine whether ear choice is associated with auditory or language brain hemispheric dominance. 
            // Can you conclude tht handedness is not independent of cell phone ear preference? Use alpha=0.05 level of significance.
            //       RE   LE  NONE
            // RH   210  103   30
            // LH    77   64   16
            std::vector<double> obs = { 210., 103., 30., 77., 64., 16. }, exp, csq;
            double n = std::accumulate(obs.begin(), obs.end(), 0.);
            // State null and alternative hypothesis: H0 = "Handedness and ear preference are independent", H1 = "Handedness and ear preference are NOT independent".
            // Find expected values using row & column totals.
            if (buildExpectedMatrix<double>(obs, exp, n, 2, 3))
            {
                // Find critical region where DoF=(rows - 1)*(cols - 1). =5.991
                unsigned DoF = 2; double cv = qchisq(1. - 0.05, DoF); print("critical value: ", cv);
                // Find chi-square. =6.744
                double cs = findChiSquare<double>(obs, exp); print("chi-square: ", cs);
                // Find p-value. =0.034
                double p = 1. - pchisq(cs, DoF); print("p-value: ", p);
                // Reject or accept H0? =reject H0
                std::cout << " We ";  DecideHypothesis(p, 0.05);
                // State conclusion in sentence.
                std::cout << " At 0.05 level of significance, there is enough evidence to conclude handedness is not independent of ear preference.\n";

            }
        }
    }


    std::cout << "Hypothesis Test for Linear Correlation (Bivariate Data)\n";
    {
        // Correlation is not causation (lurking variable). 
        {
            // Can you conclude there is a linear correlation between number of absences and the grade of a student? Use alpha=0.01 level of significance.
            // Number of absencesand the grade of a sample of students: 
            // absences = c(8, 2, 5, 12, 15, 9, 6)
            // grades = c(82, 92, 90, 51, 43, 70, 82) 
            // plot(absences, grades) 
            // lm (y, x): lm(grades~absences)
            // abline(b, m): abline(107.16, -4.213)
            // cor.test(x, y): cor.test(absences, grades)
            std::vector<double> absences = { 8, 2, 5, 12,15, 9, 6 };
            std::vector<double> grades = { 82, 92, 90, 51, 43, 70, 82 };
            std::pair<double, double> fit = lsq(absences, grades);
            std::cout << " m: " << fit.first << ", b: " << fit.second << "\n";
            double r = R(absences, grades);  print("correlation coefficeint r:", r);
            print("coefficeint of determination r^2:", pow(r, 2.));
            std::cout << " 92.4% of variation in grades (y) can be explained by variation in absences (x).\n";
            // State null and alternative hypothesis: H0 = "rho = 0", H1 = "rho != 0".
            // Find critical value, 2-tailed! =7.823
            unsigned DoF = 5 /* n - 2 */; double t = r * sqrt(DoF) / sqrt(1. - pow(r, 2.)); print("t:", t);
            // Find p-value. =0.001
            double p = (1. - pt(7.8233, DoF)) * 2.; print("p:", p);
            // Reject or accept H0? =reject H0
            std::cout << " We ";  DecideHypothesis(p, 0.02);
            // State conclusion in sentence.
            std::cout << " At 0.01 level of significance, there is enough evidence to conclude a linear correlation.\n";
        }
        {
            // Can you conclude there is a linear correlation between person's foot length and vocabulary size? Use alpha=0.05 level of significance.
            // Foot length and vocabulary size for sample of children:
            std::vector<double> foot = { 2, 5, 9, 6, 8, 8 };
            std::vector<double> vocab = { 0, 80, 124, 53, 103, 158 };
            std::pair<double, double> fit = lsq(foot, vocab);
            std::cout << " m: " << fit.first << ", b: " << fit.second << "\n";
            double r = R(vocab, foot);  print("correlation coefficeint r:", r);
            print("coefficeint of determination r^2:", pow(r, 2.));
            std::cout << " 79.9% of variation in grades (y) can be explained by variation in absences (x).\n";
            // State null and alternative hypothesis: H0 = "rho = 0", H1 = "rho != 0".
            // Find critical value, 2-tailed! = 3.991
            unsigned DoF = 4 /* n - 2 */; double t = r * sqrt(DoF) / sqrt(1. - pow(r, 2.)); print("t:", t);
            // Find p-value. =0.016
            double p = (1. - pt(3.99143, DoF)) * 2.; print("p:", p);
            // Reject or accept H0? =reject H0
            std::cout << " We ";  DecideHypothesis(p, 0.05);
            // State conclusion in sentence.
            std::cout << " At 0.05 level of significance, there is enough evidence to conclude a linear correlation.\n";
        }
    }


    // Poisson is a discrete distribution descibing a number of events in a fixed time interval/region of opportunity (0 to inf).
    std::cout << "Poisson Distributions\n";
    {
        // Number of visits to a web page follows a Poisson distribution with mean 15 visits per hour.
        // What is probability of getting 10 or less visits per hour, P(X≤10)? =0.118
        print("10 or less web visits per hour:", pPois(10, 15));
        // Equivalent to...
        //{ 
          //double p = 0.; 
            //for (int i=0; i<=10; i++) 
              //p += dPois(i, 15); 
          //print("10 or less web hits per hour:", p); 
        //} 

        // A call center gets average 13 calls per hour. What is probability that 
        // in 15 minutes call center will receive less than 6 call? =0.889
        print("Call center gets less than 6 calls:", pPois(5, 13. / 4));

        // What is probability of seeing exactly 3 blemishes on randomly selected piece of 
        // sheet metal, when on average one expects 1.2 blemishes, can be found with? =0.087
        print("3 blemishes on metal:", dPois(3, 1.2));

        // 170 oil tankers arrived randomly/independently at port over 104 days. Only 2 tankers may be unloaded per day, extra oil tankers wait.
        // What is the probability that no tankers will arrive on Tuesday? =0.195
        print("No tankers on Tuesday:", pPois(0, 170. / 104));
        // What is probability that more than 2 will arrive on Friday, so that some will wait for Saturday to be unloaded? =0.225
        print("2 tankers on Friday:", 1. - pPois(2, 170. / 104));
        // Assuming no tankers are left from Tuesday, what is probability exactly 1 tanker will be left from Wednesday and 0 will be left from Thursday? =0.142
        print("3 tankers on Wednesday:", dPois(3, 170. / 104));
    }
}

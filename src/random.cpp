#include <vector>
#include <array>
#include <random>  
#include <cmath> 
#include <algorithm> 
#include "random.h"
#include "constants.h"

// random number seed
std::random_device rd;
std::mt19937 twister(rd());
//std::mt19937 twister(1);

/**
 * Returns a random floating-point number from a uniform real distribution of 0.0 to 1.0.
 * @return The random number.
 */
double random_real() 
{
	std::uniform_real_distribution<> dist(0.0, 1.0);
	return dist(twister);
}

/**
 * Returns a random integer number from a uniform discrete distribution of a to b.
 * @param[in] a lower value of the distribution range
 * @param[in] b upper value of the distribution range
 * @return The random number.
 */
int random_discrete(int a, int b) 
{
	std::uniform_int_distribution<> dist(a, b);
	return dist(twister);
}

/**
 * Returns a random draw (non-negative integer) from the Poisson distribution with mean lambda (using normal distribution approximation when lambda > 30)
 * @param[in] lambda mean of the distribution
 * @return The random number.
 */ 
long long int random_poisson(double lambda) 
{
	long long int result;
	if (lambda < 1e-5) { 
		result = 0;
	}
	else if (lambda > 30) {
		// use normal approximation	
		std::normal_distribution<> dist(lambda, std::sqrt(lambda)); // distribution(mean, standard deviation)
		int x = std::round(dist(twister));
		result = std::max(0, x);
	}
	else {
		// sample poisson directly
		std::poisson_distribution<> dist(lambda); // distribution(mean)
		result = dist(twister);
	}

	return result;
}

/**
 * Returns a random draw (non-negative integer) from the Binomial distribution B(N,p).
 * Uses Normal and Poisson distribution approximations for large N.
 * @param[in] n number of trials
 * @param[in] p success probability
 * @return The random number.
 */ 
long long int random_binomial(long long int n, double p) 
{
	long long int result;
	if (n*p > 10 && n*(1 - p) > 10) {
		// use normal approximation
		std::normal_distribution<> dist(n*p, std::sqrt(n*p*(1 - p))); // distribution(mean, standard deviation)
		long long int x = std::round(dist(twister));
		if (x<0) x=0;
		if (x>n) x=n;
		result = x;
	}
	else if ((n > 20 && p < 0.05) || (n > 100 && n*p < 10)) {
		// use Poisson approximation
		result = random_poisson(n*p);
	}
	else if ((n > 20 && p > 0.95) || (n > 100 && n*(1-p) < 10)) {
		// use Poisson approximation
		result = n - random_poisson(n*(1 - p));
	}
	else {
		// use binomial distribution directly
		std::binomial_distribution<> dist(n, p);
		result = dist(twister);
	}

	return result;
}

/**
 * Returns a vector of outcomes from a random draw of the Multinomial distribution with N trials where each trial has a vector of probabilities probs.
 * @param[in] n 	number of trials
 * @param[in] probs vector of probabilities for each outcome
 * @return A vector of the number of successes for each outcome (in the same order as the probabilities).
 */
std::vector<long long int> random_multinomial(long long int n, const std::vector<double>& probs) 
{
	int num_outcomes = probs.size();
	double sum_p = 0.0;
	for (int i = 0; i < num_outcomes; ++i) {
		sum_p += probs[i];
	}

	long long int n_used = n;
	std::vector<long long int> result(num_outcomes, 0);
	for (int i=0; i < num_outcomes; ++i) {
		if (n_used > 0) {
			result[i] = random_binomial(n_used, probs[i] / sum_p);
			sum_p -= probs[i];
			n_used -= result[i];
		}
		else {
			result[i] = 0;
		}
	}

	return result;
}

/**
 * Returns a vector of outcomes from a random draw of the Multinomial distribution with N trials where each trial has a vector of probabilities probs.  
 * @param[in] n 	number of trials
 * @param[in] probs array of probabilities for each outcome (each genotype).
 * @return A vector of the number of successes for each outcome (in the same order as the probabilities).
 */
std::vector<long long int> random_multinomial(long long int n, const std::array<long long int, constants::num_gen>& probs) 
{
	int num_outcomes = probs.size();
	double sum_p = 0.0;
	for (int i = 0; i < num_outcomes; ++i) {
		sum_p += probs[i];
	}

	long long int n_used = n;
	std::vector<long long int> result(num_outcomes, 0);
	for (int i=0; i < num_outcomes; ++i) {
		if (n_used > 0) {
			result[i] = random_binomial(n_used, probs[i] / sum_p);
			sum_p -= probs[i];
			n_used -= result[i];
		}
		else {
			result[i] = 0;
		}
	}

	return result;
}

/**
 * Returns a vector of outcomes from a random draw of the Multinomial distribution with N trials where each trial has a vector of probabilities probs.
 * @param[in] n 	number of trials
 * @param[in] probs array of probabilities for each outcome (each age group)
 * @return A vector of the number of successes for each outcome (in the same order as the probabilities).
 */
std::vector<long long int> random_multinomial(long long int n, const std::array<double, constants::max_dev+1>& probs) 
{
	int num_outcomes = probs.size();
	double sum_p = 0.0;
	for (int i = 0; i < num_outcomes; ++i) {
		sum_p += probs[i];
	}

	long long int n_used = n;
	std::vector<long long int> result(num_outcomes, 0);
	for (int i=0; i < num_outcomes; ++i) {
		if (n_used > 0) {
			result[i] = random_binomial(n_used, probs[i] / sum_p);
			sum_p -= probs[i];
			n_used -= result[i];
		}
		else {
			result[i] = 0;
		}
	}

	return result;
}

/**
 * Returns a random draw (non-negative floating-point number) from a lognormal distribution with desired mean des_mean and desired variance des_var.
 * @param[in] des_mean 	desired mean
 * @param[in] des_var	desired variance
 * @return The random number.
 */
double random_lognormal(double des_mean, double des_var) 
{
    double mean = std::log(std::pow(des_mean, 2) / std::sqrt(std::pow(des_mean, 2) + des_var));
    double var = std::log(1.0 + (des_var / std::pow(des_mean, 2)));
    std::lognormal_distribution<> dist(mean, std::sqrt(var));
    return dist(twister);
}

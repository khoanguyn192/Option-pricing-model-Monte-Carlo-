#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
using namespace std;

// Monte Carlo pricing algorithms for Asian-style options

// Helper function to calculate cumulative normal distribution
// This function computes the cumulative distribution function for the normal distribution.
double norm_cdf(double x) {
    return 0.5 * erfc(-x * M_SQRT1_2);
}

// Function to generate linearly spaced values between start and end
// Useful for generating the time steps for the Monte Carlo simulation.
vector<double> linspace(double start, double end, int num) {
    vector<double> result;
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

/*
 * Simple Monte Carlo pricing for arithmetic Asian options
 *
 * Args:
 *     S - Initial Stock Price         sigma - Volatility
 *     K - Strike Price                CallPut - Option type specification (Call or Put)
 *     T - Time to Maturity            N - number of time steps for each MC path
 *     r - Risk-free interest rate     M - number of simulated paths
 *
 * Returns the simulated price of an arithmetic Asian-style call or put option and the associated 
 * 95% confidence interval of the estimate.
 */
// Simulates arithmetic Asian option prices using Monte Carlo and computes confidence intervals.
pair<double, pair<double, double>> arith_asian(double S, double K, double T, double r, double sigma, const string &CallPut, int M, int N) {
    vector<double> sims(M, 0.0);
    vector<double> t = linspace(0, T, N + 1);
    double dt = T / N;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(0.0, 1.0);

    for (int i = 0; i < M; ++i) {
        vector<double> W(N + 1, 0.0);
        for (int j = 1; j <= N; ++j) {
            W[j] = W[j - 1] + d(gen) * sqrt(dt);
        }

        double sum = 0.0;
        for (int j = 0; j <= N; ++j) {
            sum += S * exp((r - 0.5 * sigma * sigma) * t[j] - sigma * W[j]);
        }
        sims[i] = sum / (N + 1);
    }

    if (CallPut == "Call") {
        for (auto &sim : sims) sim = max(sim - K, 0.0);
    } else if (CallPut == "Put") {
        for (auto &sim : sims) sim = max(K - sim, 0.0);
    }

    for (auto &sim : sims) sim *= exp(-r * T);
    double mean = accumulate(sims.begin(), sims.end(), 0.0) / M;
    double variance = 0.0;
    for (const auto &sim : sims) variance += (sim - mean) * (sim - mean);
    variance /= (M - 1);
    double std_dev = sqrt(variance);
    pair<double, double> CI = {mean - 1.96 * std_dev / sqrt(M), mean + 1.96 * std_dev / sqrt(M)};

    return {mean, CI};
}

/*
 * Simple Monte Carlo pricing for geometric Asian options
 *
 * Args:
 *     S - Initial Stock Price         sigma - Volatility
 *     K - Strike Price                CallPut - Option type specification (Call or Put)
 *     T - Time to Maturity            N - number of time steps for each MC path
 *     r - Risk-free interest rate     M - number of simulated paths
 *
 * Returns the simulated price of an geometric Asian-style call or put option and the associated 
 * 95% confidence interval of the estimate.
 */
// Simulates geometric Asian option prices using Monte Carlo and computes confidence intervals.
pair<double, pair<double, double>> geom_asian(double S, double K, double T, double r, double sigma, const string &CallPut, int M, int N) {
    vector<double> sims(M, 0.0);
    vector<double> t = linspace(0, T, N + 1);
    double dt = T / N;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(0.0, 1.0);

    for (int i = 0; i < M; ++i) {
        vector<double> W(N + 1, 0.0);
        for (int j = 1; j <= N; ++j) {
            W[j] = W[j - 1] + d(gen) * sqrt(dt);
        }

        double product = 1.0;
        for (int j = 0; j <= N; ++j) {
            product *= S * exp((r - 0.5 * sigma * sigma) * t[j] - sigma * W[j]);
        }
        sims[i] = pow(product, 1.0 / (N + 1));
    }

    if (CallPut == "Call") {
        for (auto &sim : sims) sim = max(sim - K, 0.0);
    } else if (CallPut == "Put") {
        for (auto &sim : sims) sim = max(K - sim, 0.0);
    }

    for (auto &sim : sims) sim *= exp(-r * T);
    double mean = accumulate(sims.begin(), sims.end(), 0.0) / M;
    double variance = 0.0;
    for (const auto &sim : sims) variance += (sim - mean) * (sim - mean);
    variance /= (M - 1);
    double std_dev = sqrt(variance);
    pair<double, double> CI = {mean - 1.96 * std_dev / sqrt(M), mean + 1.96 * std_dev / sqrt(M)};

    return {mean, CI};
}

int main() {
    // Example usage
    // Demonstrates pricing for both arithmetic and geometric Asian options.
    double S = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;
    string CallPut = "Call";
    int M = 10000;
    int N = 100;

    auto result = arith_asian(S, K, T, r, sigma, CallPut, M, N);
    cout << "Arithmetic Asian Option Price: " << result.first << endl;
    cout << "95% Confidence Interval: [" << result.second.first << ", " << result.second.second << "]" << endl;

    auto geom_result = geom_asian(S, K, T, r, sigma, CallPut, M, N);
    cout << "Geometric Asian Option Price: " << geom_result.first << endl;
    cout << "95% Confidence Interval: [" << geom_result.second.first << ", " << geom_result.second.second << "]" << endl;

    return 0;
}

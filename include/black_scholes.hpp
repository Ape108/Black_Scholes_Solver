#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iomanip>

struct GridParams {
    double price_ceiling;
    double time_to_maturity;
    size_t num_price_steps; // determines δS
    size_t num_time_steps; // determines δT 

};  

struct MarketParams {
    double volatility;
    double risk_free_interest;
    double strike_price;
};

struct Coefficients {
    double alpha;
    double beta;
    double gamma;
};


Coefficients calculate_coeffs(
    const double& vol, 
    const double& r, 
    const double& time_to_maturity, 
    const size_t& time_steps, 
    const size_t& i
);

std::vector<double> evaluate_rhs(
    const std::vector<double>& V_known,
    const std::vector<double>& alpha,   
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    double V_bound_lower_j, double V_bound_lower_j_plus_1,
    double V_bound_upper_j, double V_bound_upper_j_plus_1);

std::vector<double> formulate_black_scholes(const GridParams& grid, const MarketParams& market);

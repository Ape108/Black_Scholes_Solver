#pragma once

#include <vector>
#include <string>
#include <map>
#include <fstream> 
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>


/// @brief Loads key-value pairs from a CSV file into a map.
/// @details Parses the CSV file assuming the first row contains column headers
///          and the second row contains corresponding values. Each header becomes
///          a key in the returned map, mapped to its value from the second row.
/// @param file Input file stream to read the CSV data from.
/// @param print If true, prints the parsed CSV data to standard output.
/// @return A map containing header-value pairs from the CSV file.
/// @throws std::runtime_error If the CSV has fewer than 2 rows or if the data row
///         has fewer columns than the header row.
std::map<std::string, std::string> data_loader(std::ifstream &file, const bool& print=false);

struct Coefficients {
    float alpha;
    float beta;
    float gamma;
};

Coefficients calculate_coeffs(const float& vol, const float& r, const float& time_to_maturity, const size_t& time_steps, const size_t& i);


std::vector<float> evaluate_rhs(
    const std::vector<float>& V_known, // Current known prices (size M + 1)
    const std::vector<float>& alpha,   // Pre-calculated coefficients
    const std::vector<float>& beta,
    const std::vector<float>& gamma,
    float V_bound_lower_j, float V_bound_lower_j_plus_1,
    float V_bound_upper_j, float V_bound_upper_j_plus_1);


struct GridParams {
    float price_ceiling; // Max price to model (theoretical price --> infinity, ceiling for computation)
    float time_to_maturity;
    size_t num_price_steps; // determines δS
    size_t num_time_steps; // determines δT 

};  

struct MarketParams {
    float volatility;
    float risk_free_interest;
    float strike_price;
};


std::vector<float> formulate_black_scholes(const GridParams& grid, const MarketParams& market);

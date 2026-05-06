#pragma once

#define _USE_MATH_DEFINES
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

/// @brief Represents the type of a financial option.
enum class OptionType
{
    Call,
    Put
};

/// @brief Contains the dimensions and boundaries for the finite difference grid.
struct GridParams {
    double price_ceiling;       ///< The maximum underlying asset price modeled on the grid.
    double time_to_maturity;    ///< Time to expiration in years (T).
    size_t num_price_steps;     ///< Number of intervals on the price axis (determines δS).
    size_t num_time_steps;      ///< Number of intervals on the time axis (determines δT).
};  

/// @brief Contains the specific market conditions and option contract details.
struct MarketParams {
    double volatility;          ///< Annualized implied volatility (σ).
    double risk_free_interest;  ///< Annualized continuous risk-free interest rate (r).
    double strike_price;        ///< The strike price of the option contract (K).
    OptionType option_type;     ///< Whether the option is a Call or a Put.
    double dividend_yield;      ///< Annualized continuous dividend yield (q).
};

/// @brief Holds the calculated Crank-Nicolson finite difference coefficients.
struct Coefficients {
    double alpha;
    double beta;
    double gamma;
};

/// @brief Calculates the alpha, beta, and gamma coefficients for a specific price step.
/// @param vol Annualized volatility.
/// @param r Risk-free interest rate.
/// @param q Dividend yield.
/// @param time_to_maturity Total time to expiration in years.
/// @param time_steps The total number of time steps (N) in the grid.
/// @param i The current price step index.
/// @return A Coefficients struct containing alpha, beta, and gamma.
Coefficients calculate_coeffs(
    double vol, 
    double r, 
    double q,
    double time_to_maturity, 
    size_t time_steps, 
    size_t i
);

/// @brief Evaluates the Right-Hand Side (RHS) of the PDE for the explicit step.
/// @param V_known The vector of option prices from the previous (known) time step.
/// @param alpha Pre-calculated alpha coefficients for the grid.
/// @param beta Pre-calculated beta coefficients for the grid.
/// @param gamma Pre-calculated gamma coefficients for the grid.
/// @param V_bound_lower_j Lower boundary condition at current time step.
/// @param V_bound_lower_j_plus_1 Lower boundary condition at previous time step.
/// @param V_bound_upper_j Upper boundary condition at current time step.
/// @param V_bound_upper_j_plus_1 Upper boundary condition at previous time step.
/// @param rhs_buffer A pre-allocated buffer where the resulting RHS values will be written.
void evaluate_rhs(
    const std::vector<double>& V_known,
    const std::vector<double>& alpha,   
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    double V_bound_lower_j, double V_bound_lower_j_plus_1,
    double V_bound_upper_j, double V_bound_upper_j_plus_1,
    std::vector<double>& rhs_buffer
);

/// @brief Validates the mathematical integrity of the grid and market parameters.
/// @param grid The GridParams struct to validate.
/// @param market The MarketParams struct to validate.
/// @throws std::invalid_argument if any parameter would cause a mathematical error.
void param_safety_check(const GridParams &grid, const MarketParams &market);

/// @brief Solves the Black-Scholes PDE using the Crank-Nicolson method.
/// @param grid The grid dimensions and boundaries.
/// @param market The market conditions and option details.
/// @return A vector representing the option's present value across all underlying price steps.
std::vector<double> formulate_black_scholes(const GridParams& grid, const MarketParams& market);

/// @brief Computes the cumulative distribution function (CDF) for the standard normal distribution.
/// @param x The standard score (z-score).
/// @return The probability that a normally distributed random variable is less than or equal to x.
double norm_cdf(double x);

/// @brief Calculates the exact price of a European option using the closed-form Black-Scholes formula.
/// @param S Current price of the underlying asset.
/// @param K Strike price.
/// @param T Time to maturity in years.
/// @param r Risk-free interest rate.
/// @param q Dividend yield.
/// @param sigma Annualized implied volatility.
/// @param type Call or Put.
/// @return The theoretical price of the European option.
double european_price(double S, double K, double T, double r, double q, double sigma, OptionType type);

/// @brief Uses Brent's Method to back-solve for Implied Volatility given a market option price.
/// @param target_price The actual market price of the option.
/// @param S Current price of the underlying asset.
/// @param K Strike price.
/// @param T Time to maturity in years.
/// @param r Risk-free interest rate.
/// @param q Dividend yield.
/// @param type Call or Put.
/// @return The implied volatility as a decimal (e.g., 0.20 for 20%).
/// @throws std::runtime_error if the root finder fails to converge.
double calculate_implied_volatility(double target_price, double S, double K, double T, double r, double q, OptionType type);
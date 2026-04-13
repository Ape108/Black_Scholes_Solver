#include "black_scholes.hpp"
#include "linear_algebra.hpp"

Coefficients calculate_coeffs(
    const double& vol, 
    const double& r, 
    const double& time_to_maturity, 
    const size_t& time_steps, 
    const size_t& i
) {
    Coefficients coeffs;
    double delta_t = time_to_maturity / time_steps;

    coeffs.alpha = (delta_t / 4.0) * (std::pow(vol, 2.0) * std::pow(i, 2.0) - r * i);
    coeffs.beta = (-delta_t / 2.0) * (std::pow(vol, 2.0) * std::pow(i, 2.0) + r);
    coeffs.gamma = (delta_t / 4.0) * (std::pow(vol, 2.0) * std::pow(i, 2.0) + r * i);

    return coeffs;
}

std::vector<double> evaluate_rhs(
    const std::vector<double>& V_known, // Current known prices (size M + 1)
    const std::vector<double>& alpha,   // Pre-calculated coefficients
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    double V_bound_lower_j, double V_bound_lower_j_plus_1,
    double V_bound_upper_j, double V_bound_upper_j_plus_1) 
{
    size_t M = V_known.size() - 1; 
    std::vector<double> rhs(M - 1, 0.0); // Matrix size is M-1

    // Loop through the internal grid nodes (i = 1 to M-1)
    for (size_t i = 1; i <= M - 1; ++i) {
        size_t rhs_idx = i - 1; // rhs is 0-indexed, grid is 1-indexed

        rhs[rhs_idx] = (1.0 + beta[i]) * V_known[i];
        
        // Add adjacent internal nodes (avoiding boundary out-of-bounds)
        if (i > 1)   rhs[rhs_idx] += alpha[i] * V_known[i - 1];
        if (i < M-1) rhs[rhs_idx] += gamma[i] * V_known[i + 1];
    }

    // Add boundaries to the first and last rows of the RHS vector
    rhs[0] += alpha[1] * (V_bound_lower_j + V_bound_lower_j_plus_1);
    rhs[M - 2] += gamma[M - 1] * (V_bound_upper_j + V_bound_upper_j_plus_1);

    return rhs;
}

std::vector<double> formulate_black_scholes(const GridParams& grid, const MarketParams& market) {
    
    if (grid.num_price_steps == 0) {
        throw std::invalid_argument("formulate_black_scholes: num_price_steps must be greater than 0");
    }
    
    if (grid.num_time_steps == 0) {
        throw std::invalid_argument("formulate_black_scholes: num_time_steps must be greater than 0");
    }
    
    if (grid.price_ceiling <= 0.0) {
        throw std::invalid_argument("formulate_black_scholes: price_ceiling must be positive");
    }
    
    if (grid.time_to_maturity <= 0.0) {
        throw std::invalid_argument("formulate_black_scholes: time_to_maturity must be positive");
    }
    
    if (market.volatility < 0.0) {
        throw std::invalid_argument("formulate_black_scholes: volatility cannot be negative");
    }
    
    if (market.strike_price <= 0.0) {
        throw std::invalid_argument("formulate_black_scholes: strike_price must be positive");
    }
    
    size_t M = grid.num_price_steps;
    size_t N = grid.num_time_steps;

    // 1. Pre-calculate coefficients for i = 0 to M
    std::vector<double> alpha(M + 1, 0.0);
    std::vector<double> beta(M + 1, 0.0);
    std::vector<double> gamma(M + 1, 0.0);

    for (size_t i = 0; i <= M; ++i) {
        Coefficients c = calculate_coeffs(market.volatility, market.risk_free_interest, grid.time_to_maturity, N, i);
        alpha[i] = c.alpha;
        beta[i] = c.beta;
        gamma[i] = c.gamma;
    }
    
    // 2. Build the left-hand Matrix A (Tridiagonal)
    std::vector<double> a_diag(M - 1); // Main diagonal: (1 - beta_i)
    std::vector<double> b_diag(M - 2); // Upper diagonal: -gamma_i
    std::vector<double> c_diag(M - 2); // Lower diagonal: -alpha_i

    for (size_t i = 1; i <= M - 1; ++i) {
        a_diag[i - 1] = 1.0 - beta[i];
        if (i < M - 1) b_diag[i - 1] = -gamma[i];
        if (i > 1)     c_diag[i - 2] = -alpha[i];
    } 

    // 3. Decompose Matrix A only once
    Decomposed LU = lu_decomposition(a_diag, b_diag, c_diag);

    // 4. Set up the terminal payoff at expiration (j = N)
    double delta_S = grid.price_ceiling / M;
    std::vector<double> V(M + 1, 0.0);
    for (size_t i = 0; i <= M; ++i) {
        double S_i = i * delta_S;
        if (market.option_type == OptionType::Call) {
            V[i] = std::max(0.0, S_i - market.strike_price);
        } else {
            V[i] = std::max(0.0, market.strike_price - S_i);
        }
    }

    // 5. Time-stepping loop (backward induction)
    double delta_t = grid.time_to_maturity / N;

    for (int j = N - 1; j >= 0; --j) {
        double t_j = j * delta_t;
        double t_j_plus_1 = (j + 1) * delta_t;

        // Calculate Boundaries for current and future step
        double V_lower_j, V_lower_j1, V_upper_j, V_upper_j1;

        if (market.option_type == OptionType::Call) {
            V_lower_j = 0.0;
            V_lower_j1 = 0.0;
            
            V_upper_j = grid.price_ceiling - market.strike_price * std::exp(-market.risk_free_interest * (grid.time_to_maturity - t_j));
            V_upper_j1 = grid.price_ceiling - market.strike_price * std::exp(-market.risk_free_interest * (grid.time_to_maturity - t_j_plus_1));
        } else { // OptionType::Put
            V_lower_j = market.strike_price * std::exp(-market.risk_free_interest * (grid.time_to_maturity - t_j));
            V_lower_j1 = market.strike_price * std::exp(-market.risk_free_interest * (grid.time_to_maturity - t_j_plus_1));
            
            V_upper_j = 0.0;
            V_upper_j1 = 0.0;
        }
        
        // Evaluate RHS
        std::vector<double> rhs = evaluate_rhs(V, alpha, beta, gamma, V_lower_j, V_lower_j1, V_upper_j, V_upper_j1);

        // Solve the system using forward and backward substitution
        std::vector<double> y = forward_substitution(LU.lower, rhs);
        std::vector<double> x = backward_substitution(LU.upper, b_diag, y);

        // Update V for the next iteration: Internal nodes become x, boundaries become deterministic equations
        for (size_t i = 1; i <= M - 1; ++i) {
            V[i] = x[i - 1];
        }
        V[0] = V_lower_j;
        V[M] = V_upper_j;
    }
    
    return V; // This vector contains the present value of the option across all price steps.
}


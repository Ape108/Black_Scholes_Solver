#include "black_scholes.hpp"
#include "linear_algebra.hpp"

Coefficients calculate_coeffs(
    double vol, 
    double r, 
    double q,
    double time_to_maturity, 
    size_t time_steps, 
    size_t i
) {
    Coefficients coeffs;
    double delta_t = time_to_maturity / time_steps;
    double drift = r - q;

    coeffs.alpha = (delta_t / 4.0) * (std::pow(vol, 2.0) * std::pow(i, 2.0) - drift * i);
    coeffs.beta = (-delta_t / 2.0) * (std::pow(vol, 2.0) * std::pow(i, 2.0) + r);
    coeffs.gamma = (delta_t / 4.0) * (std::pow(vol, 2.0) * std::pow(i, 2.0) + drift * i);

    return coeffs;
}

void evaluate_rhs(
    const std::vector<double>& V_known, // Current known prices (size M + 1)
    const std::vector<double>& alpha,   // Pre-calculated coefficients
    const std::vector<double>& beta,
    const std::vector<double>& gamma,
    double V_bound_lower_j, double V_bound_lower_j_plus_1,
    double V_bound_upper_j, double V_bound_upper_j_plus_1,
    std::vector<double>& rhs_buffer
) 
{
    size_t M = V_known.size() - 1; 

    // Loop through the internal grid nodes (i = 1 to M-1)
    for (size_t i = 1; i <= M - 1; ++i) {
        size_t rhs_idx = i - 1; // rhs is 0-indexed, grid is 1-indexed

        rhs_buffer[rhs_idx] = (1.0 + beta[i]) * V_known[i];
        
        // Add adjacent internal nodes (avoiding boundary out-of-bounds)
        if (i > 1)   rhs_buffer[rhs_idx] += alpha[i] * V_known[i - 1];
        if (i < M-1) rhs_buffer[rhs_idx] += gamma[i] * V_known[i + 1];
    }

    // Add boundaries to the first and last rows of the RHS vector
    rhs_buffer[0] += alpha[1] * (V_bound_lower_j + V_bound_lower_j_plus_1);
    rhs_buffer[M - 2] += gamma[M - 1] * (V_bound_upper_j + V_bound_upper_j_plus_1);
}

void param_safety_check(const GridParams& grid, const MarketParams& market) {
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
}

std::vector<double> formulate_black_scholes(const GridParams& grid, const MarketParams& market) {

    param_safety_check(grid, market); // Error Handling for Invalid parameters

    size_t M = grid.num_price_steps;
    size_t N = grid.num_time_steps;

    // 1. Pre-calculate coefficients for i = 0 to M
    std::vector<double> alpha(M + 1, 0.0);
    std::vector<double> beta(M + 1, 0.0);
    std::vector<double> gamma(M + 1, 0.0);

    for (size_t i = 0; i <= M; ++i) {
        Coefficients c = calculate_coeffs(market.volatility, market.risk_free_interest, market.dividend_yield, grid.time_to_maturity, N, i);
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
    std::vector<double> intrinsic_payoff(M + 1, 0.0);
    std::vector<double> V(M + 1, 0.0); // The active grid
    for (size_t i = 0; i <= M; ++i) {
        double S_i = i * delta_S;
        if (market.option_type == OptionType::Call) {
            intrinsic_payoff[i] = std::max(0.0, S_i - market.strike_price);
        } else {
            intrinsic_payoff[i] = std::max(0.0, market.strike_price - S_i);
        }

        // Seed the terminal payoff at expiration (j = N)
        V[i] = intrinsic_payoff[i];
    }

    // 5. Time-stepping loop (backward induction for American Options)
    double delta_t = grid.time_to_maturity / N;

    std::vector<double> rhs_buffer(M - 1, 0.0);
    std::vector<double> y_buffer(M - 1, 0.0);
    std::vector<double> x_buffer(M - 1, 0.0);

    for (int j = N - 1; j >= 0; --j) {
        // American Boundary Conditions (No exponential time-decay needed)
        double V_lower_j, V_lower_j1, V_upper_j, V_upper_j1;

        if (market.option_type == OptionType::Call) {
            V_lower_j = 0.0;
            V_lower_j1 = 0.0;
            V_upper_j = grid.price_ceiling - market.strike_price;
            V_upper_j1 = grid.price_ceiling - market.strike_price;
        } else { // OptionType::Put
            V_lower_j = market.strike_price;
            V_lower_j1 = market.strike_price;
            V_upper_j = 0.0;
            V_upper_j1 = 0.0;
        }
        
        evaluate_rhs(V, alpha, beta, gamma, V_lower_j, V_lower_j1, V_upper_j, V_upper_j1, rhs_buffer);
        forward_substitution(LU.lower, rhs_buffer, y_buffer);
        backward_substitution(LU.upper, b_diag, y_buffer, x_buffer);

        // ADD THE BRENNAN-SCHWARTZ CONSTRAINT
        for (size_t i = 1; i <= M - 1; ++i) {

            V[i] = std::max(x_buffer[i - 1], intrinsic_payoff[i]);
        }
        V[0] = V_lower_j;
        V[M] = V_upper_j;
    }
    
    return V; // This vector contains the present value of the option across all price steps.
}

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Closed-Form European Black-Scholes Pricer
double european_price(double S, double K, double T, double r, double q, double sigma, OptionType type) {
    if (sigma <= 0.0) return (type == OptionType::Call) ? std::max(0.0, S - K) : std::max(0.0, K - S);

    double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    if (type == OptionType::Call) {
        return S * std::exp(-q * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
    } else {
        return K * std::exp(-r * T) * norm_cdf(-d2) - S * std::exp(-q * T) * norm_cdf(-d1);
    }
}

double calculate_implied_volatility(double target_price, double S, double K, double T, double r, double q, OptionType type) {
    // 1. Define the bracket [a, b]
    double a = 1e-4; // Lower bound (0.01%)
    double b = 5.0;  // Upper bound (500%)
    
    // Function to evaluate the error at a given volatility
    auto eval_error = [&](double vol) {
        return european_price(S, K, T, r, q, vol, type) - target_price;
    };

    double fa = eval_error(a);
    double fb = eval_error(b);

    // If the target price is completely outside our 500% volatility bounds
    if (fa * fb > 0.0) {
        throw std::runtime_error("Implied volatility root is not bracketed. Market price may be invalid or arbitrageable.");
    }

    // Brent's Method Setup
    double c = a, fc = fa;
    double d = b - a, e = d;
    double tol = 1e-6; // Precision tolerance
    int max_iter = 100;

    for (int iter = 0; iter < max_iter; ++iter) {
        if (std::abs(fc) < std::abs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double tol1 = 2.0 * 2.2204460492503131e-16 * std::abs(b) + 0.5 * tol; // Machine epsilon scaling
        double xm = 0.5 * (c - b);

        if (std::abs(xm) <= tol1 || fb == 0.0) {
            return b; // Root found!
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            double s = fb / fa;
            double p, q_val;
            
            if (a == c) {
                // Secant Method
                p = 2.0 * xm * s;
                q_val = 1.0 - s;
            } else {
                // Inverse Quadratic Interpolation
                q_val = fa / fc;
                double r_val = fb / fc;
                p = s * (2.0 * xm * q_val * (q_val - r_val) - (b - a) * (r_val - 1.0));
                q_val = (q_val - 1.0) * (r_val - 1.0) * (s - 1.0);
            }

            if (p > 0.0) q_val = -q_val;
            p = std::abs(p);

            double min1 = 3.0 * xm * q_val - std::abs(tol1 * q_val);
            double min2 = std::abs(e * q_val);

            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q_val; // Accept interpolation
            } else {
                d = xm; e = d; // Fallback to Bisection
            }
        } else {
            d = xm; e = d; // Fallback to Bisection
        }

        a = b; fa = fb;
        if (std::abs(d) > tol1) b += d;
        else b += (xm > 0.0 ? tol1 : -tol1);

        fb = eval_error(b);
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa;
            e = d = b - a;
        }
    }

    throw std::runtime_error("Brent's Method failed to converge within maximum iterations.");
}

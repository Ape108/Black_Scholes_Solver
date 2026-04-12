#include "black_scholes.hpp"
#include "linear_algebra.hpp"


std::map<std::string, std::string> data_loader(std::ifstream &file, bool print) {
    
    if (!file.is_open()) {
        throw std::runtime_error("CSV data_loader: input file stream is not open");
    }
    
    std::string line;
    std::vector<std::vector<std::string>> rows;

    while(std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }
        
        if (!row.empty()) {  // Skip empty lines
            rows.push_back(row);
        }
    }

    if (rows.size() < 2) {
        throw std::runtime_error("CSV data_loader: file must contain at least a header row and one data row, got " + std::to_string(rows.size()) + " row(s)");
    }

    const auto &header = rows[0];
    const auto &values = rows[1];
    
    if (header.empty()) {
        throw std::runtime_error("CSV data_loader: header row is empty");
    }
    
    if (values.size() < header.size()) {
        throw std::runtime_error("CSV data_loader: data row has " + std::to_string(values.size()) + " columns but header has " + std::to_string(header.size()) + " columns");
    }

    if (print) {
        std::cout << std::endl;
        // Print the parsed data
        for (size_t i = 0; i < header.size(); ++i) {
            std::cout << header[i] << ": " << values[i] << std::endl;
        }
        std::cout << std::endl;
    }

    std::map<std::string, std::string> params;
    for (size_t i=0; i<header.size(); i++) {
        params[header[i]] = values[i];
    }

    return params;
}

std::map<std::string, std::string> open_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening file:" + filename);
    }

    // Load option data from file into config
    std::map<std::string, std::string> params = data_loader(file, true);

    file.close();
    return params;
};



double price_option(
    std::vector<double>& V, 
    const size_t& num_price_steps,
    const std::map<std::string, std::string>& params
) {

    if (V.empty()) {
        throw std::invalid_argument("price_option: solution vector V cannot be empty");
    }
    
    if (num_price_steps == 0) {
        throw std::invalid_argument("price_option: num_price_steps must be greater than 0");
    }
    
    if (V.size() != num_price_steps + 1) {
        throw std::invalid_argument("price_option: solution vector size (" + std::to_string(V.size()) + 
                                   ") must equal num_price_steps + 1 (" + std::to_string(num_price_steps + 1) + ")");
    }

    double current_price = std::stod(params.at("current_price"));
    // Extracting the specific option price
    double price_ceiling = std::stod(params.at("strike_price")) * 2.0;
    
    if (price_ceiling <= 0.0) {
        throw std::invalid_argument("price_option: price ceiling must be positive");
    }
    
    double delta_S = price_ceiling / num_price_steps;

    // Find where our current price lands on the grid
    double exact_idx = current_price / delta_S;

    // Get the integer indices immediately below and above it
    size_t lower_idx = static_cast<size_t>(std::floor(exact_idx));
    size_t upper_idx = lower_idx + 1;
    
    // Validate array bounds
    if (upper_idx >= V.size()) {
        throw std::runtime_error("price_option: calculated index " + std::to_string(upper_idx) + 
                                " exceeds solution vector size " + std::to_string(V.size()));
    }

    // Calculate how close the price is to the upper node (weighting)
    double weight_upper = exact_idx - lower_idx;
    double weight_lower = 1.0 - weight_upper;

    // Interpolate the exact option price
    double option_price = (V[lower_idx] * weight_lower) + (V[upper_idx] * weight_upper);
    
    return option_price;
}

void print_option_diff(const double& option_price, const std::map<std::string, std::string>& params) {
    double actual_price = std::stod(params.at("actual_price"));
    double underlier_price = std::stod(params.at("current_price"));

    std::cout << "\n=========================================" << std::endl;
    std::cout << " Underlier Price:   $" << underlier_price << std::endl;
    std::cout << " Calculated Option: $" << std::fixed << std::setprecision(4) << option_price << std::endl;
    std::cout << " Actual Option: $" << std::fixed << std::setprecision(4) << actual_price << std::endl;
    std::cout << "=========================================\n" << std::endl;
}

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
        V[i] = std::max(0.0, S_i - market.strike_price); // Assuming Call Option
    }

    // 5. Time-stepping loop (backward induction)
    double delta_t = grid.time_to_maturity / N;

    for (int j = N - 1; j >= 0; --j) {
        double t_j = j * delta_t;
        double t_j_plus_1 = (j + 1) * delta_t;

        // Calculate Call Option boundaries for current and future step
        double V_lower_j = 0.0;
        double V_lower_j1 = 0.0;

        double V_upper_j = grid.price_ceiling - market.strike_price * std::exp(-market.risk_free_interest * (grid.time_to_maturity - t_j));
        double V_upper_j1 = grid.price_ceiling - market.strike_price * std::exp(-market.risk_free_interest * (grid.time_to_maturity - t_j_plus_1));

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

std::vector<double> evaluate_system(
    const size_t& num_time_steps, 
    const size_t& num_price_steps, 
    const std::map<std::string, std::string>& params
) {
    
    if (num_time_steps == 0) {
        throw std::invalid_argument("evaluate_system: num_time_steps must be greater than 0");
    }
    
    if (num_price_steps == 0) {
        throw std::invalid_argument("evaluate_system: num_price_steps must be greater than 0");
    }
    
    GridParams grid; 
    MarketParams market;

    grid.num_price_steps = num_price_steps;
    grid.num_time_steps = num_time_steps;
    
    try {
        grid.price_ceiling = std::stod(params.at("strike_price")) * 2.0;
        grid.time_to_maturity = std::stod(params.at("T"));
        market.risk_free_interest = std::stod(params.at("risk_free_rate"));
        market.strike_price = std::stod(params.at("strike_price"));
        market.volatility = std::stod(params.at("implied_vol"));
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("evaluate_system: required parameter missing from CSV: " + std::string(e.what()));
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("evaluate_system: invalid numeric value in CSV: " + std::string(e.what()));
    }

    std::vector<double> V = formulate_black_scholes(grid, market);

    return V;
}
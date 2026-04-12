#include "black_scholes.hpp"
#include "linear_algebra.hpp"

int main() {

    std::string filename = "option_chain.csv";
    
    size_t NUM_PRICESTEPS = 760;
    size_t NUM_TIMESTEPS = 200;
    
    // Validate input parameters
    if (NUM_PRICESTEPS == 0 || NUM_TIMESTEPS == 0) {
        std::cerr << "Error: NUM_PRICESTEPS and NUM_TIMESTEPS must be greater than 0" << std::endl;
        return 1;
    }

    try {
        // Load and parse CSV file
        std::map<std::string, std::string> params = open_file(filename);
        
        // Solve the Black-Scholes PDE using finite difference method
        std::vector<double> solution_vector = evaluate_system(NUM_TIMESTEPS, NUM_PRICESTEPS, params);
        
        // Calculate the option price at the current spot price
        double option_price = price_option(solution_vector, NUM_PRICESTEPS, params);
        
        // Display results
        print_option_diff(option_price, params);
        
        return 0;
    } 
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument error: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::runtime_error& e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range error: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        return 1;
    }
}
#include "black_scholes.hpp"
#include "linear_algebra.hpp"

int main() {

    std::string filename = "/workspaces/Black_Scholes_Solver/option_chain.csv";
    
    size_t NUM_PRICESTEPS = 760;
    size_t NUM_TIMESTEPS = 200;

    try {
        std::map<std::string, std::string> params = open_file(filename);
        std::vector<float> solution_vector = evaluate_system(NUM_TIMESTEPS, NUM_PRICESTEPS, params);
        // Run Black-Scholes Pricing Calculation
        float option_price = price_option(solution_vector, NUM_PRICESTEPS, params);
        print_option_diff(option_price, params);

    } 
    catch (const std::runtime_error& e) {
        std::cerr << "Caught: " << e.what() << std::endl;
    }
    return 0;
}
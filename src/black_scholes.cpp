#include "black_scholes.hpp"
#include <stdexcept>

std::map<std::string, std::string> data_loader(std::ifstream &file, const bool& print) {
    
    std::string line;
    std::vector<std::vector<std::string>> rows;

    while(std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }
        rows.push_back(row);
    }

    if (rows.size() < 2) {
        throw std::runtime_error("CSV data_loader: file must contain at least header and one data row");
    }

    const auto &header = rows[0];
    const auto &values = rows[1];
    if (values.size() < header.size()) {
        throw std::runtime_error("CSV data_loader: data row has fewer columns than header");
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
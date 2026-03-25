#pragma once

#include <vector>
#include <string>
#include <map>
#include <fstream> 
#include <sstream>
#include <iostream>

// I'll write a data loader to grab the inputs from the CSV 
// Then I will discretize the grid,
// and initialize the boundary conditions.

std::map<std::string, std::string> data_loader(std::ifstream &file);
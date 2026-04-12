# Black-Scholes PDE Solver

A high-performance Python library for pricing financial options. The core engine is written in C++ and approximates the Black-Scholes Partial-Differential Equation using the Crank-Nicolson Finite Difference Method. The resulting tridiagonal linear system is solved in O(N) time complexity using the Thomas Algorithm.

This project leverages `pybind11` and `scikit-build-core` to expose the C++ numerical engine directly to Python, allowing for rapid data acquisition (e.g., via `yfinance`) alongside microsecond-level PDE solving.

## Installation

To build and install the library directly from the source code, ensure you have a modern C++ compiler and CMake installed, navigate to the root directory, and run:

```bash
pip install .
```

## Quick Start

Once installed, you can import the compiled C++ engine directly into any Python script.

```python
import black_scholes_solver

# 1. Define the grid parameters
grid = black_scholes_solver.GridParams()
grid.price_ceiling = 380.0
grid.time_to_maturity = 0.0082
grid.num_price_steps = 760
grid.num_time_steps = 200

# 2. Define the market parameters
market = black_scholes_solver.MarketParams()
market.volatility = 0.2754
market.risk_free_interest = 0.0359
market.strike_price = 190.0

# 3. Solve the PDE (Executes natively in C++)
# Returns a Python list containing the option value at every price step
V = black_scholes_solver.formulate_black_scholes(grid, market)
```

## References

SkanCity Academy. (2023, October 6). 🟢05 - Thomas Algorithm for Solving Tri-diagonal Matrix Systems [Video]. YouTube. https://www.youtube.com/watch?v=vzqwV-REmkw

Smolski, A. (2023, December 30). Crank-Nicholson (Finite Difference) with Black-Scholes (with code). Medium. https://antonismolski.medium.com/crank-nicholson-with-black-scholes-with-code-a27c0df17555

Zientziateka. (2019, May 31). Matrix representation of the Crank-Nicholson method for the Black-Scholes equation [Video]. YouTube. https://youtu.be/5mp-2zqo6hY?si=rUIu-qep44Q8UE0L

Zientziateka. (2019, May 31). The Crank-Nicholson method for the Black-Scholes equation [Video]. YouTube. https://youtu.be/XHa81xxpj6I?si=ftj2lmpp0gn7GNto
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
market.option_type = black_scholes_solver.OptionType.Call

# 3. Solve the PDE (Executes natively in C++)
# Returns a Python list containing the option value at every price step
V = black_scholes_solver.formulate_black_scholes(grid, market)
```

### Calculating Implied Volatility
The library also includes a high-speed root-finder using Brent's Method to calculate arbitrage-free implied volatility directly from market prices.

```python
import black_scholes_solver

iv = black_scholes_solver.calculate_implied_volatility(
    target_price=2.5450,
    S=196.50,
    K=197.50,
    T=0.0081,
    r=0.0360,
    q=0.0,
    type=black_scholes_solver.OptionType.Call
)
print(f"Implied Volatility: {iv:.4f}")
```

![volatility_smile.png]

## Performance

The numerical engine is built to handle computationally intensive stochastic grid calculations with minimal latency. By bypassing the Python Global Interpreter Lock (GIL) and executing raw 64-bit C++ machine code, the library achieves a roughly **30x execution speedup** over native Python implementations.

**Benchmark (2000 Price Steps x 2000 Time Steps):**
* **Pure Python:** ~2.51 seconds
* **Compiled C++:** ~0.08 seconds

## Limitations & Mathematical Assumptions

While this engine is built for microsecond execution and high-precision PDE solving, it currently relies on several standard quantitative assumptions that may introduce minor drift when compared to institutional pricing feeds (e.g., OptionMetrics, Bloomberg):

* **European IV Approximation for American Options:** The `calculate_implied_volatility` root-finder utilizes Brent's Method over the closed-form European Black-Scholes equation. Because American options contain an early-exercise premium, using a European formula to extract IV from an American market price will slightly artificially inflate the resulting volatility. (Note: The PDE *does* properly price American options via the Brennan-Schwartz constraint, but the fast root-finder assumes European exercise).
* **Continuous Dividend Yields:** The solver models dividends as a continuous annualized yield ($q$). It does not currently support discrete dividend schedules (lumpy cash flows on specific ex-dividend dates). For underlyings with massive, irregular dividends, this continuous approximation will cause slight pricing drift.
* **Constant Interest Rates:** The `MarketParams` struct accepts a single, scalar constant for the risk-free rate ($r$). The engine does not natively support a full yield curve or term structure of interest rates. 
* **Mid-Price Illiquidity:** The implied volatility pipeline is optimized to target the exact Bid-Ask Mid-Price. In highly illiquid options with blown-out spreads, the arithmetic midpoint may not represent the true market clearing price, which can cause Brent's Method to map an exaggerated volatility smirk.


## References

Brennan, M. J., & Schwartz, E. S. (1977). The valuation of American put options. The Journal of Finance, 32(2), 449–462. https://www.jstor.org/stable/2326779

SkanCity Academy. (2023, October 6). 🟢05 - Thomas Algorithm for Solving Tri-diagonal Matrix Systems [Video]. YouTube. https://www.youtube.com/watch?v=vzqwV-REmkw

Smolski, A. (2023, December 30). Crank-Nicholson (Finite Difference) with Black-Scholes (with code). Medium. https://antonismolski.medium.com/crank-nicholson-with-black-scholes-with-code-a27c0df17555

Zientziateka. (2019, May 31). Matrix representation of the Crank-Nicholson method for the Black-Scholes equation [Video]. YouTube. https://youtu.be/5mp-2zqo6hY?si=rUIu-qep44Q8UE0L

Zientziateka. (2019, May 31). The Crank-Nicholson method for the Black-Scholes equation [Video]. YouTube. https://youtu.be/XHa81xxpj6I?si=ftj2lmpp0gn7GNto
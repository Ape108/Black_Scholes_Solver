import math
import time

import black_scholes_solver

# 1. PURE PYTHON IMPLEMENTATION 

def lu_decomposition(a, b, c):
    n = len(a)
    lower = []
    upper = [a[0]]

    for i in range(1, n):
        denominator = upper[i-1]
        l_i = c[i-1] / denominator
        lower.append(l_i)
        u_i = a[i] - (l_i * b[i-1])
        upper.append(u_i)
    return lower, upper

def forward_substitution(lower, b):
    y = [b[0]]
    n = len(lower)
    for i in range(1, n + 1):
        y_i = b[i] - (lower[i-1] * y[i-1])
        y.append(y_i)
    return y

def backward_substitution(u, b, y):
    n = len(u)
    stack = [y[n-1] / u[n-1]]
    
    for i in range(n-1, 0, -1):
        x_i = y[i-1] - (b[i-1] * stack[-1])
        x_i /= u[i-1]
        stack.append(x_i)
        
    stack.reverse()
    return stack

# Added 'q' (dividend yield) and generalized drift
def calculate_coeffs(vol, r, q, time_to_maturity, time_steps, i):
    delta_t = time_to_maturity / time_steps
    drift = r - q
    
    alpha = (delta_t / 4.0) * ((vol**2) * (i**2) - drift * i)
    beta = (-delta_t / 2.0) * ((vol**2) * (i**2) + r)
    gamma = (delta_t / 4.0) * ((vol**2) * (i**2) + drift * i)
    return alpha, beta, gamma

def evaluate_rhs(V_known, alpha, beta, gamma, V_lower_j, V_lower_j1, V_upper_j, V_upper_j1):
    M = len(V_known) - 1
    rhs = [0.0] * (M - 1)

    for i in range(1, M):
        rhs_idx = i - 1
        rhs[rhs_idx] = (1.0 + beta[i]) * V_known[i]
        
        if i > 1:
            rhs[rhs_idx] += alpha[i] * V_known[i - 1]
        if i < M - 1:
            rhs[rhs_idx] += gamma[i] * V_known[i + 1]

    rhs[0] += alpha[1] * (V_lower_j + V_lower_j1)
    rhs[M - 2] += gamma[M - 1] * (V_upper_j + V_upper_j1)
    return rhs

# Added 'q' and 'option_type'
def naive_formulate_black_scholes(price_ceiling, time_to_maturity, num_price_steps, num_time_steps, vol, r, q, strike, option_type):
    M = num_price_steps
    N = num_time_steps

    alpha = [0.0] * (M + 1)
    beta = [0.0] * (M + 1)
    gamma = [0.0] * (M + 1)

    for i in range(M + 1):
        a, b, g = calculate_coeffs(vol, r, q, time_to_maturity, N, i)
        alpha[i] = a
        beta[i] = b
        gamma[i] = g

    a_diag = [0.0] * (M - 1)
    b_diag = [0.0] * (M - 2)
    c_diag = [0.0] * (M - 2)

    for i in range(1, M):
        a_diag[i - 1] = 1.0 - beta[i]
        if i < M - 1: b_diag[i - 1] = -gamma[i]
        if i > 1:     c_diag[i - 2] = -alpha[i]

    lower, upper = lu_decomposition(a_diag, b_diag, c_diag)

    delta_S = price_ceiling / M
    V = [0.0] * (M + 1)
    
    # Terminal Payoff handles Calls and Puts
    for i in range(M + 1):
        S_i = i * delta_S
        if option_type == 'Call':
            V[i] = max(0.0, S_i - strike)
        else:
            V[i] = max(0.0, strike - S_i)

    # Time-stepping loop (backward induction)
    for j in range(N - 1, -1, -1):
        
        # American Boundary Conditions
        if option_type == 'Call':
            V_lower_j = 0.0
            V_lower_j1 = 0.0
            V_upper_j = price_ceiling - strike
            V_upper_j1 = price_ceiling - strike
        else:
            V_lower_j = strike
            V_lower_j1 = strike
            V_upper_j = 0.0
            V_upper_j1 = 0.0

        rhs = evaluate_rhs(V, alpha, beta, gamma, V_lower_j, V_lower_j1, V_upper_j, V_upper_j1)
        y = forward_substitution(lower, rhs)
        x = backward_substitution(upper, b_diag, y)

        # Brennan-Schwartz American Early Exercise Constraint
        for i in range(1, M):
            S_i = i * delta_S
            if option_type == 'Call':
                intrinsic_value = max(0.0, S_i - strike)
            else:
                intrinsic_value = max(0.0, strike - S_i)
                
            V[i] = max(x[i - 1], intrinsic_value)
            
        V[0] = V_lower_j
        V[M] = V_upper_j

    return V


# 2. BENCHMARK RACE

def run_benchmark():
    # Set expensive grid parameters to stress-test both engines
    ITERATIONS = 5
    PRICE_CEILING = 400.0
    TTM = 1.5
    PRICE_STEPS = 2000
    TIME_STEPS = 2000
    VOL = 0.25
    RATE = 0.05
    DIV_YIELD = 0.015 # Added dividend yield parameter
    STRIKE = 200.0

    print(f"Benchmarking Crank-Nicolson American PDE Solver...")
    print(f"Grid Size: {PRICE_STEPS} Price Steps x {TIME_STEPS} Time Steps")
    print(f"Running {ITERATIONS} iterations to calculate average latency...\n")

    # Python Test
    print("Running Pure Python Implementation...")
    py_start = time.perf_counter()
    for _ in range(ITERATIONS):
        # Passing DIV_YIELD and 'Call' as option type
        py_V = naive_formulate_black_scholes(
            PRICE_CEILING, TTM, PRICE_STEPS, TIME_STEPS, VOL, RATE, DIV_YIELD, STRIKE, 'Call'
        )
    py_total = time.perf_counter() - py_start
    py_avg = py_total / ITERATIONS
    print(f"Python Average Execution Time: {py_avg:.5f} seconds\n")

    # Compiled C++ Test
    print("Running Compiled C++ Library...")
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = PRICE_CEILING
    grid.time_to_maturity = TTM
    grid.num_price_steps = PRICE_STEPS
    grid.num_time_steps = TIME_STEPS

    market = black_scholes_solver.MarketParams()
    market.volatility = VOL
    market.risk_free_interest = RATE
    market.dividend_yield = DIV_YIELD  # Updating the new struct property
    market.strike_price = STRIKE
    market.option_type = black_scholes_solver.OptionType.Call # Native Enum

    cpp_start = time.perf_counter()
    for _ in range(ITERATIONS):
        cpp_V = black_scholes_solver.formulate_black_scholes(grid, market)
    cpp_total = time.perf_counter() - cpp_start
    cpp_avg = cpp_total / ITERATIONS
    print(f"C++ Average Execution Time:    {cpp_avg:.5f} seconds\n")

    # Results & Sanity Check
    diff = abs(py_V[PRICE_STEPS // 2] - cpp_V[PRICE_STEPS // 2])

    speedup = py_avg / cpp_avg
    print("-" * 40)
    print(f"RESULTS:")
    print(f"Math Drift: ${diff:.8f}")
    print(f"Speedup:    {speedup:.1f}x FASTER than Python")
    print("-" * 40)

if __name__ == "__main__":
    run_benchmark()
import black_scholes_solver

def test_call_option_execution():
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = 380.0
    grid.time_to_maturity = 0.0082
    grid.num_price_steps = 100
    grid.num_time_steps = 50

    market = black_scholes_solver.MarketParams()
    market.volatility = 0.2754
    market.risk_free_interest = 0.0359
    market.strike_price = 190.0
    market.option_type = black_scholes_solver.OptionType.Call # ADDED

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    assert len(V) == grid.num_price_steps + 1
    assert V[-1] > 0.0 # Call option has value at the ceiling
    assert V[0] == 0.0 # Call option is worthless at stock price of $0

def test_put_option_execution():
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = 380.0
    grid.time_to_maturity = 0.0082
    grid.num_price_steps = 100
    grid.num_time_steps = 50

    market = black_scholes_solver.MarketParams()
    market.volatility = 0.2754
    market.risk_free_interest = 0.0359
    market.strike_price = 190.0
    market.option_type = black_scholes_solver.OptionType.Put # NEW PUT TEST

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    assert len(V) == grid.num_price_steps + 1
    assert V[-1] == 0.0 # Put option is worthless at the ceiling (infinity)
    assert V[0] > 0.0   # Put option has maximum value when stock price is $0

def test_implied_volatility_calculation():
    # Known market conditions
    S = 100.0
    K = 100.0
    T = 1.0
    r = 0.05
    q = 0.02
    target_volatility = 0.25 # We expect it to find 25%
    
    # The mathematical European price for these exact parameters is roughly $10.027
    target_price = 10.0271 
    
    calculated_iv = black_scholes_solver.calculate_implied_volatility(
        target_price=target_price,
        S=S,
        K=K,
        T=T,
        r=r,
        q=q,
        type=black_scholes_solver.OptionType.Call
    )
    
    # Assert that Brent's method found the 25% volatility to within 4 decimal places
    assert abs(calculated_iv - target_volatility) < 1e-4

def test_implied_volatility_put_option():
    # Test that it works for puts as well
    S = 50.0
    K = 55.0
    T = 0.5
    r = 0.04
    q = 0.0
    
    target_price = 5.8559 # Expected price for a Put at 30% vol
    
    calculated_iv = black_scholes_solver.calculate_implied_volatility(
        target_price=target_price,
        S=S,
        K=K,
        T=T,
        r=r,
        q=q,
        type=black_scholes_solver.OptionType.Put
    )
    
    assert abs(calculated_iv - 0.30) < 1e-4
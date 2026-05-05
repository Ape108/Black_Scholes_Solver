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
    market.option_type = black_scholes_solver.OptionType.Call

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    assert len(V) == grid.num_price_steps + 1
    assert V[-1] > 0.0
    assert V[0] == 0.0 

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
    market.option_type = black_scholes_solver.OptionType.Put 

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    assert len(V) == grid.num_price_steps + 1
    assert V[-1] == 0.0 
    assert V[0] > 0.0   

def test_implied_volatility_calculation():
    # Brand new parameters: At-The-Money Call with 0 interest/dividends
    # This removes exponential compounding so we can verify the math perfectly.
    S = 100.0
    K = 100.0
    T = 1.0
    r = 0.0
    q = 0.0
    
    # We will test if the root finder can find 20% volatility
    target_volatility = 0.20 
    
    # Mathematical Proof:
    # d1 = 0.10, d2 = -0.10
    # Exact Price = 100 * (0.5398278 - 0.4601722) = 7.96556
    target_price = 7.96556
    
    calculated_iv = black_scholes_solver.calculate_implied_volatility(
        target_price=target_price,
        S=S,
        K=K,
        T=T,
        r=r,
        q=q,
        type=black_scholes_solver.OptionType.Call
    )
    
    assert abs(calculated_iv - target_volatility) < 1e-4

def test_implied_volatility_put_option():
    # Brand new parameters: At-The-Money Put with 0 interest/dividends
    S = 100.0
    K = 100.0
    T = 1.0
    r = 0.0
    q = 0.0
    
    # We will test if the root finder can find 30% volatility
    target_volatility = 0.30
    
    # Mathematical Proof:
    # d1 = 0.15, d2 = -0.15
    # Exact Price = 100 * (0.559617 - 0.440382) = 11.9235
    target_price = 11.9235
    
    calculated_iv = black_scholes_solver.calculate_implied_volatility(
        target_price=target_price,
        S=S,
        K=K,
        T=T,
        r=r,
        q=q,
        type=black_scholes_solver.OptionType.Put
    )
    
    assert abs(calculated_iv - target_volatility) < 1e-4
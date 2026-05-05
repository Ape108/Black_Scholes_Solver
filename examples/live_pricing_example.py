import yfinance as yf
import pandas as pd
from datetime import datetime, timezone
import black_scholes_solver

def get_risk_free_rate() -> float:
    """Fetches the current risk-free rate from the 13-week Treasury Bill (^IRX)."""
    irx = yf.Ticker("^IRX")
    return round((irx.history(period='1d')['Close'].iloc[-1] / 100.0), 4)

def calculate_time_to_maturity(target_expiry: str) -> float:
    """Calculates the time to maturity (T) in years, standardized to 365 days."""
    expiry_date = datetime.strptime(target_expiry, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    today = datetime.now(timezone.utc)
    
    # Options typically expire at 4:00 PM EST, roughly 20:00 UTC
    expiry_exact = expiry_date.replace(hour=20, minute=0, second=0)
    fractional_days = (expiry_exact - today).total_seconds() / 86400.0
    
    return round(max(fractional_days / 365.0, 1e-5), 6)

def interpolate_price(V: list, current_price: float, price_ceiling: float, num_price_steps: int) -> float:
    """Interpolates the exact option price from the PDE grid."""
    delta_S = price_ceiling / num_price_steps
    exact_idx = current_price / delta_S
    lower_idx = int(exact_idx)
    upper_idx = lower_idx + 1

    weight_upper = exact_idx - lower_idx
    weight_lower = 1.0 - weight_upper

    return (V[lower_idx] * weight_lower) + (V[upper_idx] * weight_upper)

def main():
    target_stock = "NVDA"
    print(f"\nTicker: {target_stock}", flush=True)

    # --- 1. Data Acquisition ---
    ticker = yf.Ticker(target_stock)
    current_price = round(ticker.history(period='1d')['Close'].iloc[-1], 2)
    risk_free_rate = get_risk_free_rate()
    div_yield = ticker.info.get('dividendYield', 0.0) or 0.0

    # Get Option Chain (Index 1 to avoid 0DTE math errors)
    target_expiry = ticker.options[1] 
    chain = ticker.option_chain(target_expiry)
    
    # Find the closest At-The-Money (ATM) Call Option
    calls = chain.calls
    target_call = calls.iloc[(calls['strike'] - current_price).abs().argsort()[:1]].iloc[0]

    # --- 2. Parameter Extraction ---
    strike_price = round(target_call['strike'], 2)
    yfinance_iv = round(target_call['impliedVolatility'], 4)
    bid = round(target_call['bid'], 4)
    ask = round(target_call['ask'], 4)
    actual_price = round((bid + ask) / 2.0, 4)
    T = calculate_time_to_maturity(target_expiry)

    print(f"Market Bid: ${bid:.2f} | Market Ask: ${ask:.2f}")
    print(f"True Market Mid-Price: ${actual_price:.4f}\n")

    # --- 3. Implied Volatility Calculation ---
    calculated_iv = black_scholes_solver.calculate_implied_volatility(
        target_price=actual_price,
        S=current_price,
        K=strike_price,
        T=T,
        r=risk_free_rate,
        q=div_yield,
        type=black_scholes_solver.OptionType.Call
    )

    print(f"yFinance Quoted IV:   {yfinance_iv:.4f}")
    print(f"C++ Calculated IV:    {calculated_iv:.4f}\n")

    print(f"--- Extracted Parameters ---")
    print(f"S0 (Price): ${current_price:.2f} | K (Strike): ${strike_price} | T (Years): {T:.4f}")
    # Note: Updated print statement to reflect the calculated_iv, not the stale yFinance IV
    print(f"Vol (sigma): {calculated_iv:.4f} | Risk-Free Rate (r): {risk_free_rate:.4f}\n")

    # --- 4. PDE Grid Execution ---
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = strike_price * 2.0
    grid.time_to_maturity = T
    grid.num_price_steps = 2500
    grid.num_time_steps = 2500

    market = black_scholes_solver.MarketParams()
    market.volatility = calculated_iv
    market.risk_free_interest = risk_free_rate
    market.strike_price = strike_price
    market.option_type = black_scholes_solver.OptionType.Call
    market.dividend_yield = div_yield

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    # --- 5. Final Output ---
    theoretical_price = interpolate_price(V, current_price, grid.price_ceiling, grid.num_price_steps)
    print(f"Theoretical C++ Price: ${theoretical_price:.4f}")

if __name__ == "__main__":
    main()
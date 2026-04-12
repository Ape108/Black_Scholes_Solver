import yfinance as yf
import pandas as pd
from datetime import datetime, timezone
import bs_solver # My Compiled C++ Library 

def main():
    start_time = datetime.now()

    target_stock = "NVDA"

    ticker = yf.Ticker(target_stock)
    print("\nTicker: " + target_stock, flush=True)

    # Current stock price (S0)
    current_price = round(ticker.history(period='1d')['Close'].iloc[-1], 2)

    # Risk-free interest rate
    irx = yf.Ticker("^IRX")
    risk_free_rate = round((irx.history(period='1d')['Close'].iloc[-1] / 100.0), 4)

    expiration_dates = ticker.options
    target_expiry = expiration_dates[1] # Using index 1 to avoid 0DTE math errors (T=0).
    chain = ticker.option_chain(target_expiry)

    calls = chain.calls
    target_call = calls.iloc[(calls['strike'] - current_price).abs().argsort()[:1]].iloc[0]

    strike_price = round(target_call['strike'], 2)
    implied_vol = round(target_call['impliedVolatility'], 4)

    # Grab the bid and ask instead of the lastPrice
    bid = round(target_call['bid'], 4)
    ask = round(target_call['ask'], 4)

    # Calculate the true market consensus (Mid-Price)
    actual_price = round((bid + ask) / 2.0, 4)

    print(f"Market Bid: ${bid:.2f} | Market Ask: ${ask:.2f}")
    print(f"True Market Mid-Price: ${actual_price:.4f}\n")

    # Calculate time to maturity (T) in years
    expiry_date = datetime.strptime(target_expiry, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    today = datetime.now(timezone.utc)
    days_to_expiry = (expiry_date - today).days
    T = round(max(days_to_expiry / 365.0, 1e-5), 4) # Prevent division by zero if it expires today

    print(f"--- Extracted Parameters ---")
    print(f"S0 (Price): ${current_price:.2f} | K (Strike): ${strike_price} | T (Years): {T:.4f}")
    print(f"Vol (sigma): {implied_vol:.4f} | Risk-Free Rate (r): {risk_free_rate:.4f}\n")

    end_time = datetime.now()

    elapsed_time = end_time - start_time
    print("Data extraction took", elapsed_time, "seconds.\n")

    start_time = datetime.now()
    # 1. Instantiate the C++ structs directly in Python
    grid = bs_solver.GridParams()
    grid.price_ceiling = strike_price * 2.0
    grid.time_to_maturity = T
    grid.num_price_steps = 760
    grid.num_time_steps = 200

    market = bs_solver.MarketParams()
    market.volatility = implied_vol
    market.risk_free_interest = risk_free_rate
    market.strike_price = strike_price

    # 2. Call the C++ engine (Executes in compiled C++ speed)
    # pybind11 automatically converts the returned std::vector<float> into a Python list
    V = bs_solver.formulate_black_scholes(grid, market)

    # 3. Interpolate the exact price using Python
    delta_S = grid.price_ceiling / grid.num_price_steps
    exact_idx = current_price / delta_S
    lower_idx = int(exact_idx)
    upper_idx = lower_idx + 1

    weight_upper = exact_idx - lower_idx
    weight_lower = 1.0 - weight_upper

    theoretical_price = (V[lower_idx] * weight_lower) + (V[upper_idx] * weight_upper)

    print(f"Theoretical C++ Price: ${theoretical_price:.4f}")

    end_time = datetime.now()
    elapsed_time = end_time - start_time
    print("\nCalculation took", elapsed_time, "seconds.\n")

if __name__ == "__main__":
    main()
import yfinance as yf
import pandas as pd
from datetime import datetime, timezone

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

# Calculate time to maturity (T) in years
expiry_date = datetime.strptime(target_expiry, '%Y-%m-%d').replace(tzinfo=timezone.utc)
today = datetime.now(timezone.utc)
days_to_expiry = (expiry_date - today).days
T = round(max(days_to_expiry / 365.0, 1e-5), 4) # Prevent division by zero if it expires today

print(f"--- Extracted Parameters ---")
print(f"S0 (Price): ${current_price:.2f} | K (Strike): ${strike_price} | T (Years): {T:.4f}")
print(f"Vol (sigma): {implied_vol:.4f} | Risk-Free Rate (r): {risk_free_rate:.4f}\n")

df = pd.DataFrame(data={"ticker":[target_stock], "current_price":[current_price], "strike_price":[strike_price],
                             "T":[T], "implied_vol":[implied_vol], "risk_free_rate":[risk_free_rate]})

df.to_csv('option_chain.csv', index=False)
# %% Imports
import yfinance as yf

# %% Sample
target_stock = "NVDA"

dat = yf.Ticker(target_stock)

print(dat.info)
print(dat.calendar)
print(dat.analyst_price_targets)
print(dat.quarterly_income_stmt)
print(dat.history(period='1mo'))
print(dat.option_chain(dat.options[0]).calls)


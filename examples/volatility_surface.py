import yfinance as yf
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timezone
import black_scholes_solver

def get_risk_free_rate() -> float:
    irx = yf.Ticker("^IRX")
    return round((irx.history(period='1d')['Close'].iloc[-1] / 100.0), 4)

def calculate_time_to_maturity(target_expiry: str) -> float:
    expiry_date = datetime.strptime(target_expiry, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    today = datetime.now(timezone.utc)
    expiry_exact = expiry_date.replace(hour=20, minute=0, second=0)
    fractional_days = (expiry_exact - today).total_seconds() / 86400.0
    return round(max(fractional_days / 365.0, 1e-5), 6)

def plot_volatility_smile(df: pd.DataFrame, current_price: float, target_stock: str):
    """Renders a volatility smile chart comparing API data to C++ calculated data."""
    # Ensure matplotlib uses a clean, professional style
    plt.style.use('dark_background') # Or 'seaborn-v0_8-darkgrid' for a light theme
    plt.figure(figsize=(12, 7))

    # Plot the yFinance benchmark (Dashed Gray)
    plt.plot(df['Strike'], df['yFinance_IV'] * 100, 
             label='yFinance API IV', 
             color='#888888', 
             linestyle='--', 
             linewidth=2)

    # Plot your C++ Calculated IV (Solid Cyan/Blue)
    plt.plot(df['Strike'], df['CPP_IV'] * 100, 
             label='C++ Calculated IV', 
             color='#00e5ff', 
             linewidth=2.5)

    # Add a vertical line to show where the stock is currently trading (At-The-Money)
    plt.axvline(x=current_price, 
                color='#ff3366', 
                linestyle=':', 
                linewidth=2, 
                label=f'Current Price (ATM): ${current_price:.2f}')

    # Formatting the chart
    plt.title(f"{target_stock} Volatility Skew / Smile", fontsize=16, pad=15)
    plt.xlabel("Strike Price ($)", fontsize=12)
    plt.ylabel("Implied Volatility (%)", fontsize=12)
    
    plt.grid(True, alpha=0.2)
    plt.legend(fontsize=12)
    plt.tight_layout()

    # Render the window
    plt.savefig("volatility_smile.png")

def main():
    target_stock = "NVDA"
    print(f"\nFetching options chain for: {target_stock}...")

    ticker = yf.Ticker(target_stock)
    current_price = round(ticker.history(period='1d')['Close'].iloc[-1], 2)
    risk_free_rate = get_risk_free_rate()
    div_yield = ticker.info.get('dividendYield', 0.0) or 0.0

    target_expiry = ticker.options[1] 
    chain = ticker.option_chain(target_expiry)
    calls = chain.calls
    T = calculate_time_to_maturity(target_expiry)

    print(f"S0: ${current_price:.2f} | T: {T:.4f} yrs | r: {risk_free_rate:.4f}\n")
    print(f"{'Strike':<10} | {'Mid-Price':<10} | {'yFinance IV':<15} | {'C++ IV':<15}")
    print("-" * 60)

    results = []

    # Iterate through every strike in the call chain
    for index, row in calls.iterrows():
        strike = float(row['strike'])
        yfinance_iv = float(row['impliedVolatility'])
        
        bid = float(row['bid'])
        ask = float(row['ask'])
        
        # Skip severely illiquid options where bid/ask spread is broken
        if bid <= 0.0 or ask <= 0.0:
            continue
            
        mid_price = (bid + ask) / 2.0

        try:
            # Blast the parameters into your C++ engine
            calc_iv = black_scholes_solver.calculate_implied_volatility(
                target_price=mid_price,
                S=current_price,
                K=strike,
                T=T,
                r=risk_free_rate,
                q=div_yield,
                type=black_scholes_solver.OptionType.Call
            )
            
            print(f"${strike:<9.2f} | ${mid_price:<9.2f} | {yfinance_iv:<15.4f} | {calc_iv:<15.4f}")
            
            results.append({
                "Strike": strike,
                "MidPrice": mid_price,
                "yFinance_IV": yfinance_iv,
                "CPP_IV": calc_iv
            })
            
        except Exception as e:
            # Brent's method will fail if the mid_price implies an arbitrage opportunity
            # (e.g., option is priced lower than its intrinsic value)
            pass

    # You can now convert 'results' to a pandas DataFrame to plot your volatility smile
    df_surface = pd.DataFrame(results)

    # Check if we successfully calculated any rows before plotting
    if not df_surface.empty:
        print("\nRendering Volatility Smile Chart...")
        plot_volatility_smile(df_surface, current_price, target_stock)
    else:
        print("\nNo valid options data to plot.")
    
if __name__ == "__main__":
    main()
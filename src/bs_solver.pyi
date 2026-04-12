
class GridParams:
    price_ceiling: float
    time_to_maturity: float
    num_price_steps: int
    num_time_steps: int
    def __init__(self) -> None: ...

class MarketParams:
    volatility: float
    risk_free_interest: float
    strike_price: float
    def __init__(self) -> None: ...

def formulate_black_scholes(grid: GridParams, market: MarketParams) -> list[float]: ...
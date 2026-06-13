# Brownian Motion

**A Rust library for simulating Brownian motion and Geometric Brownian Motion (GBM)** — the mathematical foundation of options pricing, risk modeling, and stochastic process simulation.

## Why It Matters

Brownian motion (the Wiener process) is the random walk that underlies much of quantitative finance. The Black-Scholes options pricing model assumes stock prices follow Geometric Brownian Motion: `dS = μS dt + σS dW`, where `dW` is a Brownian increment.

Key financial applications:
- **Monte Carlo pricing** — simulate thousands of price paths to estimate exotic option values
- **Value at Risk (VaR)** — estimate portfolio loss distributions
- **Stress testing** — simulate worst-case price trajectories
- **Algorithmic trading** — backtest strategies against synthetic price data

Brownian motion also appears in physics (particle diffusion), biology (molecular motors), and engineering (signal noise).

The **quadratic variation** of a Brownian path equals its time length — this is the mathematical property that makes stochastic calculus (Itô's lemma) work. The **maximum drawdown** metric measures the worst peak-to-trough decline in a simulated path — a key risk indicator.

## How It Works

**Brownian increments**: A single step of Brownian motion is `√dt · Z` where Z is a standard normal random variable. The library generates Z using the **Box-Muller transform**: `Z = √(−2·ln(U₁)) · cos(2π·U₂)` where U₁, U₂ are uniform [0,1) random variables from the `rand` crate.

**Path generation**: Each `path(n)` call generates n+1 points starting at 0.0, adding an increment at each step. The cumulative sum forms a random walk.

**Geometric Brownian Motion**: Models an asset price as `S(t+dt) = S(t) · exp((μ − σ²/2)·dt + σ·dW)`. The drift term `μ − σ²/2` accounts for the lognormal distribution of prices (Itô's correction). The volatility σ scales the random component.

**Quadratic variation**: Computed as `Σ(Δx)²` over the path. For a true Brownian path, this equals the time horizon — a fundamental property of the Wiener process.

**Maximum drawdown**: Tracks the running peak and computes the largest percentage decline from any peak to a subsequent trough.

## Quick Start

```rust
use brownian_motion::{BrownianMotion, GeometricBrownianMotion, quadratic_variation, max_drawdown};

// Standard Brownian motion (Wiener process)
let bm = BrownianMotion::new(0.01); // dt = 0.01
let path = bm.path(1000); // 1000 steps
println!("Final position: {:.4}", path.last().unwrap());
println!("Quadratic variation: {:.4}", quadratic_variation(&path));

// Geometric Brownian Motion — stock price simulation
let gbm = GeometricBrownianMotion::new(0.05, 0.2, 1.0/252.0); // μ=5%, σ=20%, daily
let prices = gbm.path(100.0, 252); // $100 start, 252 trading days
let risk = max_drawdown(&prices);
println!("Max drawdown: {:.1}%", risk * 100.0);
```

## API

- **`BrownianMotion`** — Standard Wiener process: `new(dt)`, `increment()`, `path(n)`
- **`GeometricBrownianMotion`** — GBM for asset prices: `new(mu, sigma, dt)`, `path(s0, n)`
- **`quadratic_variation(path)` → `f64`** — Sum of squared increments
- **`max_drawdown(prices)` → `f64`** — Largest peak-to-trough decline (0.0–1.0)

## Architecture Notes

Provides the stochastic process simulation primitives for SuperInstance quantitative analysis tools. Used in portfolio risk assessment, options pricing, and synthetic data generation. See the [architecture overview](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT

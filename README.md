# brownian-motion

A Rust library for **stochastic process simulation**, implementing standard Brownian motion (Wiener process), Geometric Brownian Motion (GBM), and path analytics including quadratic variation and maximum drawdown.

## Why It Matters

Brownian motion is the mathematical foundation of quantitative finance, physical diffusion, and stochastic optimization. The Wiener process W(t) underlies:

- **Black–Scholes option pricing** — GBM models asset price dynamics
- **Kalman filtering** — process noise is Brownian
- **Molecular dynamics** — Langevin equations add friction to Brownian paths
- **Stochastic gradient descent** — gradient noise is modeled as Brownian

Monte Carlo simulation of these processes requires efficient, correct incremental sampling — exactly what this crate provides.

## How It Works

### Standard Wiener Process

A Wiener process W(t) satisfies three axioms:

1. W(0) = 0 (almost surely)
2. W(t) − W(s) ~ N(0, t − s) for t > s
3. Independent increments

The incremental simulation uses the **Box–Muller transform** to generate standard normals:

$$Z = \sqrt{-2 \ln U_1} \cdot \cos(2\pi U_2)$$

Each increment is scaled: $\Delta W = Z \cdot \sqrt{\Delta t}$

This ensures the variance scales correctly: $\text{Var}(\Delta W) = \Delta t$.

### Geometric Brownian Motion

GBM follows the Itô SDE:

$$dS_t = \mu S_t \, dt + \sigma S_t \, dW_t$$

The Euler–Maruyama discretization (used here) is:

$$S_{t+\Delta t} = S_t \cdot \exp\left(\left(\mu - \frac{\sigma^2}{2}\right)\Delta t + \sigma \Delta W\right)$$

The drift correction $-\frac{\sigma^2}{2}$ arises from Itô's lemma — this is the difference between Itô and Stratonovich calculus.

### Quadratic Variation

For a sampled path, the quadratic variation is:

$$[W]_T = \sum_{i=1}^{n} (W_{t_i} - W_{t_{i-1}})^2$$

For a true Brownian path, $[W]_T = T$ almost surely as the mesh $\max \Delta t_i \to 0$. This is the key distinguishing feature of Brownian motion (vs. differentiable functions, where QV = 0).

### Maximum Drawdown

$$\text{MDD} = \max_t \frac{\text{peak}(t) - S(t)}{\text{peak}(t)}$$

computed in O(n) with a single pass tracking the running peak.

### Big-O Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `path(n)` | O(n) | O(n) |
| `increment()` | O(1) | O(1) |
| `quadratic_variation(&path)` | O(n) | O(1) |
| `max_drawdown(&prices)` | O(n) | O(1) |

## Quick Start

```rust
use brownian_motion::{BrownianMotion, GeometricBrownianMotion, quadratic_variation, max_drawdown};

// Standard Wiener process
let bm = BrownianMotion::new(0.01); // dt = 0.01
let path = bm.path(1000);           // 1000 steps
assert_eq!(path[0], 0.0);

// Geometric Brownian Motion (Black-Scholes dynamics)
let gbm = GeometricBrownianMotion::new(0.05, 0.2, 1.0/252.0); // μ=5%, σ=20%, daily
let prices = gbm.path(100.0, 252);  // 1 year of daily prices

// Analytics
let qv = quadratic_variation(&path);
let mdd = max_drawdown(&prices);
```

## API

| Type | Method | Description |
|------|--------|-------------|
| `BrownianMotion` | `new(dt: f64)` | Create with time step |
| `BrownianMotion` | `increment() → f64` | Single ΔW sample |
| `BrownianMotion` | `path(n: usize) → Vec<f64>` | Full path of n steps |
| `GeometricBrownianMotion` | `new(mu, sigma, dt)` | Create with drift & volatility |
| `GeometricBrownianMotion` | `path(s0, n) → Vec<f64>` | Price path from S₀ |
| `quadratic_variation(&[f64]) → f64` | Sum of squared increments |
| `max_drawdown(&[f64]) → f64` | Peak-to-trough decline fraction |

## Architecture Notes

The **γ + η = C** link: the Box–Muller transform (γ) generates Gaussian increments from uniform randomness, while the √dt scaling (η) calibrates them to the correct variance. Together they conserve the Wiener invariant C — the resulting process has the correct covariance structure Cov(W(s), W(t)) = min(s, t), verified by the quadratic variation converging to T.

## References

- Wiener, N. (1923). *Differential-space.* Journal of Mathematics and Physics, 2(1–4), 131–174.
- Itô, K. (1944). *Stochastic Integral.* Proc. Imperial Acad. Tokyo, 20, 519–524.
- Black, F., & Scholes, M. (1973). *The Pricing of Options and Corporate Liabilities.* JPE, 81(3), 637–654.
- Kloeden, P. E., & Platen, E. (1992). *Numerical Solution of Stochastic Differential Equations.* Springer.
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering.* Springer.

## License

MIT

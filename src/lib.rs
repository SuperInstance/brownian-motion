//! Brownian Motion Simulation
//!
//! Geometric Brownian motion and standard Wiener process utilities.

use std::f64::consts::PI;

/// Standard Wiener process (Brownian motion) simulator.
pub struct BrownianMotion {
    dt: f64,
}

impl BrownianMotion {
    pub fn new(dt: f64) -> Self {
        assert!(dt > 0.0, "time step must be positive");
        Self { dt }
    }

    /// Generate a single Brownian increment (scaled normal).
    pub fn increment(&self) -> f64 {
        let u1: f64 = rand::random();
        let u2: f64 = rand::random();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        z * self.dt.sqrt()
    }

    /// Simulate a full path of n steps, returning cumulative positions.
    pub fn path(&self, n: usize) -> Vec<f64> {
        let mut positions = Vec::with_capacity(n + 1);
        positions.push(0.0);
        let mut x = 0.0;
        for _ in 0..n {
            x += self.increment();
            positions.push(x);
        }
        positions
    }
}

/// Geometric Brownian Motion with drift μ and volatility σ.
pub struct GeometricBrownianMotion {
    mu: f64,
    sigma: f64,
    dt: f64,
}

impl GeometricBrownianMotion {
    pub fn new(mu: f64, sigma: f64, dt: f64) -> Self {
        assert!(sigma > 0.0, "volatility must be positive");
        assert!(dt > 0.0, "time step must be positive");
        Self { mu, sigma, dt }
    }

    /// Simulate a GBM path starting from s0 for n steps.
    pub fn path(&self, s0: f64, n: usize) -> Vec<f64> {
        let bm = BrownianMotion::new(self.dt);
        let mut prices = Vec::with_capacity(n + 1);
        prices.push(s0);
        let mut s = s0;
        for _ in 0..n {
            let dw = bm.increment();
            s *= (self.mu - 0.5 * self.sigma * self.sigma) * self.dt + self.sigma * dw;
            s = s.max(0.0);
            prices.push(s);
        }
        prices
    }
}

/// Compute the quadratic variation of a path.
pub fn quadratic_variation(path: &[f64]) -> f64 {
    path.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum()
}

/// Compute the maximum drawdown of a price path.
pub fn max_drawdown(prices: &[f64]) -> f64 {
    let mut peak = 0.0_f64;
    let mut max_dd = 0.0_f64;
    for &p in prices {
        if p > peak {
            peak = p;
        }
        let dd = (peak - p) / peak.max(1e-10);
        if dd > max_dd {
            max_dd = dd;
        }
    }
    max_dd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_brownian_path_starts_at_zero() {
        let bm = BrownianMotion::new(0.01);
        let p = bm.path(100);
        assert_eq!(p[0], 0.0);
        assert_eq!(p.len(), 101);
    }

    #[test]
    fn test_gbm_path_length() {
        let gbm = GeometricBrownianMotion::new(0.05, 0.2, 1.0 / 252.0);
        let p = gbm.path(100.0, 252);
        assert_eq!(p.len(), 253);
    }

    #[test]
    fn test_quadratic_variation() {
        let path = vec![0.0, 1.0, 2.0, 3.0];
        assert_eq!(quadratic_variation(&path), 3.0);
    }
}

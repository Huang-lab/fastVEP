//! SVM evaluation for de novo donor splice site prediction.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Pre-loaded SVM model for de novo donor prediction.
pub struct SvmModel {
    /// Feature names in order.
    pub feature_names: Vec<String>,
    /// Support vectors: each row is [feature_values..., alpha].
    pub support_vectors: Vec<Vec<f64>>,
    /// Feature centering values.
    pub center: HashMap<String, f64>,
    /// Feature scaling values.
    pub scale: HashMap<String, f64>,
    /// SVM parameter: decision boundary offset.
    pub rho: f64,
    /// Platt scaling parameter A.
    pub prob_a: f64,
    /// Platt scaling parameter B.
    pub prob_b: f64,
    /// RBF kernel gamma parameter.
    pub gamma: f64,
}

impl SvmModel {
    /// Load SVM model from a directory containing sv.txt, center.txt, scale.txt, misc.txt.
    pub fn load(dir: &Path) -> Result<Self, String> {
        // Load support vectors
        let sv_path = dir.join("sv.txt");
        let file = std::fs::File::open(&sv_path)
            .map_err(|e| format!("Cannot open {}: {}", sv_path.display(), e))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // First line is header with feature names
        let header_line = lines
            .next()
            .ok_or("Empty sv.txt")?
            .map_err(|e| format!("Read error: {}", e))?;
        let feature_names: Vec<String> = header_line
            .split('\t')
            .filter(|s| !s.is_empty())
            .map(|s| s.trim().to_string())
            .collect();
        // Last column is "alpha", feature names are all but last
        let n_features = feature_names.len() - 1; // exclude alpha

        let mut support_vectors = Vec::new();
        for line in lines {
            let line = line.map_err(|e| format!("Read error: {}", e))?;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let values: Vec<f64> = trimmed
                .split('\t')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if values.len() == n_features + 1 {
                support_vectors.push(values);
            }
        }

        // Load center, scale, misc
        let center = load_kv_file(&dir.join("center.txt"))?;
        let scale = load_kv_file(&dir.join("scale.txt"))?;
        let misc = load_kv_file(&dir.join("misc.txt"))?;

        let rho = *misc.get("rho").ok_or("Missing rho in misc.txt")?;
        let prob_a = *misc.get("probA").ok_or("Missing probA in misc.txt")?;
        let prob_b = *misc.get("probB").ok_or("Missing probB in misc.txt")?;
        let gamma = *misc.get("gamma").ok_or("Missing gamma in misc.txt")?;

        Ok(Self {
            feature_names: feature_names[..n_features].to_vec(),
            support_vectors,
            center,
            scale,
            rho,
            prob_a,
            prob_b,
            gamma,
        })
    }

    /// Evaluate the SVM decision function for the given features.
    /// Returns the raw margin value.
    pub fn decision_function(&self, features: &HashMap<String, f64>) -> f64 {
        // Scale features
        let scaled: Vec<f64> = self
            .feature_names
            .iter()
            .map(|name| {
                let val = features.get(name).copied().unwrap_or(0.0);
                let center = self.center.get(name).copied().unwrap_or(0.0);
                let scale = self.scale.get(name).copied().unwrap_or(1.0);
                (val - center) / scale
            })
            .collect();

        let n_features = self.feature_names.len();
        let mut margin = 0.0;

        for sv in &self.support_vectors {
            let sv_features = &sv[..n_features];
            let alpha = sv[n_features];

            // RBF kernel
            let k = rbf_kernel(&scaled, sv_features, self.gamma);
            margin += alpha * k;
        }

        margin - self.rho
    }

    /// Convert SVM margin to probability using Platt scaling.
    pub fn predict_probability(&self, margin: f64) -> f64 {
        let logit = self.prob_a * margin + self.prob_b;
        1.0 / (1.0 + (-logit).exp())
    }

    /// Evaluate SVM and return probability.
    pub fn predict(&self, features: &HashMap<String, f64>) -> f64 {
        let margin = self.decision_function(features);
        self.predict_probability(margin)
    }
}

/// RBF (Radial Basis Function) kernel.
fn rbf_kernel(v1: &[f64], v2: &[f64], gamma: f64) -> f64 {
    let sq_dist: f64 = v1
        .iter()
        .zip(v2.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    (-gamma * sq_dist).exp()
}

fn load_kv_file(path: &Path) -> Result<HashMap<String, f64>, String> {
    let file = std::fs::File::open(path)
        .map_err(|e| format!("Cannot open {}: {}", path.display(), e))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() == 2 {
            if let Ok(val) = parts[1].parse() {
                map.insert(parts[0].to_string(), val);
            }
        }
    }
    Ok(map)
}

/// Logistic regression scoring (used for disruption probability).
pub fn logreg(features: &HashMap<String, f64>, coefficients: &HashMap<String, f64>) -> f64 {
    let mut logit = coefficients.get("(Intercept)").copied().unwrap_or(0.0);
    for (name, &coef) in coefficients {
        if name == "(Intercept)" {
            continue;
        }
        let val = features.get(name).copied().unwrap_or(0.0);
        logit += coef * val;
    }
    1.0 / (1.0 + (-logit).exp())
}

/// Load logistic regression coefficients from a file (name\tvalue per line).
pub fn load_logreg_model(path: &Path) -> Result<HashMap<String, f64>, String> {
    load_kv_file(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rbf_kernel_same() {
        let v = vec![1.0, 2.0, 3.0];
        let k = rbf_kernel(&v, &v, 0.25);
        assert!((k - 1.0).abs() < 1e-10); // same vector → distance=0, exp(0)=1
    }

    #[test]
    fn test_rbf_kernel_different() {
        let v1 = vec![0.0, 0.0];
        let v2 = vec![1.0, 1.0];
        let k = rbf_kernel(&v1, &v2, 0.25);
        // sq_dist = 2.0, exp(-0.25*2.0) = exp(-0.5) ≈ 0.6065
        assert!((k - 0.6065).abs() < 0.001);
    }

    #[test]
    fn test_logreg() {
        let mut features = HashMap::new();
        features.insert("x".to_string(), 1.0);
        let mut coefficients = HashMap::new();
        coefficients.insert("(Intercept)".to_string(), 0.0);
        coefficients.insert("x".to_string(), 0.0);
        let prob = logreg(&features, &coefficients);
        assert!((prob - 0.5).abs() < 1e-10); // logit=0 → prob=0.5
    }

    #[test]
    fn test_logreg_strong_positive() {
        let mut features = HashMap::new();
        features.insert("x".to_string(), 10.0);
        let mut coefficients = HashMap::new();
        coefficients.insert("(Intercept)".to_string(), 0.0);
        coefficients.insert("x".to_string(), 1.0);
        let prob = logreg(&features, &coefficients);
        assert!(prob > 0.999); // logit=10 → very high probability
    }
}

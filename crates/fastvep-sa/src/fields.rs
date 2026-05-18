//! Field type definitions for the .osa2 format.
//!
//! Each annotation field (e.g., gnomAD AF, ClinVar significance) has a type
//! that determines how it's stored as a u32 integer and reconstructed for output.

use serde::{Deserialize, Serialize};

/// How a field's values are stored internally.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FieldType {
    /// Raw integer (u32). Optional zigzag for signed values.
    Integer,
    /// Float scaled to integer: `stored = (value * multiplier) as u32`.
    Float,
    /// String mapped to u32 index via a lookup table.
    Categorical,
    /// Boolean flag: 0 or 1.
    Flag,
    /// Opaque JSON string stored separately (not in u32 arrays).
    JsonBlob,
}

/// Definition of a single annotation field within an .osa2 file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Field {
    /// Original VCF INFO field name (e.g., "AF").
    pub field: String,
    /// Output alias (e.g., "gnomad_af"). Used as the JSON key.
    pub alias: String,
    /// How the field is stored.
    pub ftype: FieldType,
    /// For Float fields: multiply by this before casting to u32.
    /// E.g., 2_000_000 for allele frequencies (6 significant digits).
    #[serde(default = "default_multiplier")]
    pub multiplier: u32,
    /// Whether to apply zigzag encoding (for signed values).
    #[serde(default)]
    pub zigzag: bool,
    /// Sentinel value stored when the field is missing. Default: u32::MAX.
    #[serde(default = "default_missing")]
    pub missing_value: u32,
    /// For Categorical fields: the string to use for missing values.
    #[serde(default = "default_missing_string")]
    pub missing_string: String,
    /// Human-readable description.
    #[serde(default)]
    pub description: String,
}

fn default_multiplier() -> u32 { 1 }
fn default_missing() -> u32 { u32::MAX }
fn default_missing_string() -> String { ".".into() }

impl Field {
    /// Encode a float value as a u32 using this field's multiplier.
    ///
    /// Non-finite inputs (NaN, +/-Inf) — and scalings that overflow to
    /// non-finite — map to the field's `missing_value`. Finite scaled values
    /// are clamped into the storable range: `[0, u32::MAX]` for non-zigzag
    /// (or `[0, u32::MAX - 1]` when the missing sentinel is `u32::MAX`) and
    /// `[i32::MIN, i32::MAX]` for zigzag (with `i32::MIN` excluded when the
    /// missing sentinel is `u32::MAX`, since its zigzag encoding collides).
    /// Without this guard, `f64 as u32` would silently produce 0 for NaN and
    /// saturate to `u32::MAX` for huge values — which collides with the
    /// default missing sentinel and corrupts data.
    #[inline]
    pub fn encode_float(&self, value: f64) -> u32 {
        if !value.is_finite() {
            return self.missing_value;
        }
        let scaled = value * self.multiplier as f64;
        if !scaled.is_finite() {
            return self.missing_value;
        }
        if self.zigzag {
            // Clamp into i32 range before zigzag-encoding signed values.
            // When the missing sentinel is u32::MAX, reserve it by excluding
            // i32::MIN, whose zigzag encoding is also u32::MAX.
            let min_storable = if self.missing_value == u32::MAX {
                (i32::MIN + 1) as f64
            } else {
                i32::MIN as f64
            };
            let clamped = scaled.clamp(min_storable, i32::MAX as f64) as i32;
            crate::zigzag::encode(clamped)
        } else {
            // Saturate in u32 range. Avoid colliding with the missing sentinel
            // when the field uses u32::MAX as its missing value.
            let max_storable = if self.missing_value == u32::MAX {
                (u32::MAX - 1) as f64
            } else {
                u32::MAX as f64
            };
            let bounded = scaled.clamp(0.0, max_storable);
            bounded as u32
        }
    }

    /// Decode a stored u32 back to a float.
    #[inline]
    pub fn decode_float(&self, stored: u32) -> f64 {
        if stored == self.missing_value {
            return f64::NAN;
        }
        let raw = if self.zigzag {
            crate::zigzag::decode(stored) as f64
        } else {
            stored as f64
        };
        let m = self.multiplier as f64;
        if m == 0.0 { raw } else { raw / m }
    }

    /// Encode an integer value, optionally with zigzag.
    ///
    /// For zigzag fields, the value is clamped into `[i32::MIN, i32::MAX]`
    /// before encoding (with `i32::MIN` excluded when the missing sentinel is
    /// `u32::MAX`, since its zigzag encoding collides with the sentinel).
    ///
    /// For non-zigzag fields, negative inputs or values exceeding `u32::MAX`
    /// are returned as `missing_value` rather than silently wrapping into a
    /// large unrelated `u32`.
    #[inline]
    pub fn encode_int(&self, value: i64) -> u32 {
        if self.zigzag {
            // When the missing sentinel is u32::MAX, zigzag-encoding i32::MIN
            // would collide with that sentinel. Clamp away from i32::MIN so we
            // never need to remap the encoded value onto another valid extreme.
            let min = if self.missing_value == u32::MAX {
                (i32::MIN as i64) + 1
            } else {
                i32::MIN as i64
            };
            let clamped = value.clamp(min, i32::MAX as i64) as i32;
            crate::zigzag::encode(clamped)
        } else {
            // Treat negatives or overflow as missing rather than silently wrapping.
            if value < 0 || value > u32::MAX as i64 {
                self.missing_value
            } else {
                value as u32
            }
        }
    }

    /// Decode a stored u32 back to an integer.
    #[inline]
    pub fn decode_int(&self, stored: u32) -> i64 {
        if stored == self.missing_value {
            return 0;
        }
        if self.zigzag {
            crate::zigzag::decode(stored) as i64
        } else {
            stored as i64
        }
    }
}

/// Format a decoded value as a JSON fragment string.
pub fn format_value(field: &Field, stored: u32, strings: Option<&[String]>) -> String {
    if stored == field.missing_value {
        return "null".into();
    }
    match field.ftype {
        FieldType::Float => {
            let val = field.decode_float(stored);
            if val.abs() < 1e-10 { "0".into() }
            else { format!("{:.6e}", val) }
        }
        FieldType::Integer => {
            let val = field.decode_int(stored);
            format!("{}", val)
        }
        FieldType::Categorical => {
            if let Some(strs) = strings {
                let idx = stored as usize;
                if idx < strs.len() {
                    format!("\"{}\"", strs[idx])
                } else {
                    "null".into()
                }
            } else {
                format!("{}", stored)
            }
        }
        FieldType::Flag => {
            if stored != 0 { "true".into() } else { "false".into() }
        }
        FieldType::JsonBlob => {
            // JsonBlob values are stored separately, not in u32 arrays
            "null".into()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_float_round_trip() {
        let field = Field {
            field: "AF".into(), alias: "af".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let original = 0.001234;
        let stored = field.encode_float(original);
        let decoded = field.decode_float(stored);
        assert!((decoded - original).abs() < 1e-6);
    }

    #[test]
    fn test_zigzag_float() {
        let field = Field {
            field: "score".into(), alias: "score".into(), ftype: FieldType::Float,
            multiplier: 10_000, zigzag: true, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let original = -2.5;
        let stored = field.encode_float(original);
        let decoded = field.decode_float(stored);
        assert!((decoded - original).abs() < 0.001);
    }

    #[test]
    fn test_missing_value() {
        let field = Field {
            field: "AF".into(), alias: "af".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        assert!(field.decode_float(u32::MAX).is_nan());
        assert_eq!(format_value(&field, u32::MAX, None), "null");
    }

    #[test]
    fn test_encode_nan_and_inf_become_missing() {
        let field = Field {
            field: "AF".into(), alias: "af".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        assert_eq!(field.encode_float(f64::NAN), u32::MAX);
        assert_eq!(field.encode_float(f64::INFINITY), u32::MAX);
        assert_eq!(field.encode_float(f64::NEG_INFINITY), u32::MAX);
    }

    #[test]
    fn test_encode_saturation_does_not_collide_with_missing() {
        // Without bounds, value * multiplier may saturate to u32::MAX and
        // collide with the default missing sentinel.
        let field = Field {
            field: "x".into(), alias: "x".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let encoded = field.encode_float(1e30);
        assert_ne!(encoded, u32::MAX, "saturation must not collide with missing");
        assert!(encoded < u32::MAX);
    }

    #[test]
    fn test_encode_int_negative_without_zigzag_treated_missing() {
        let field = Field {
            field: "n".into(), alias: "n".into(), ftype: FieldType::Integer,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        assert_eq!(field.encode_int(-5), u32::MAX);
        assert_eq!(field.encode_int(0), 0);
        assert_eq!(field.encode_int(42), 42);
    }

    #[test]
    fn test_categorical() {
        let field = Field {
            field: "sig".into(), alias: "sig".into(), ftype: FieldType::Categorical,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: String::new(),
        };
        let strings = vec!["Benign".to_string(), "Pathogenic".to_string()];
        assert_eq!(format_value(&field, 0, Some(&strings)), "\"Benign\"");
        assert_eq!(format_value(&field, 1, Some(&strings)), "\"Pathogenic\"");
    }
}

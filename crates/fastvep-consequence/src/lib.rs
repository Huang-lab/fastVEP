mod edits;
mod predictor;
mod splice;
pub mod sv_predictor;

pub use predictor::{
    AlleleConsequenceResult, ConsequencePredictor, Flank, PredictionResult, ProteinWindow,
    TranscriptConsequence,
};

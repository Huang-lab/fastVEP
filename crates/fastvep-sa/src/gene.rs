//! Gene-level annotation reader/writer (.oga files).
//!
//! Used for gene-keyed databases (OMIM, gnomAD gene scores, ClinGen).
//! Annotations are looked up by gene symbol.

use crate::common::{GeneRecord, MAX_INDEX_PAYLOAD, OGA_MAGIC, SCHEMA_VERSION};
use anyhow::Result;
use fastvep_cache::annotation::GeneAnnotationProvider;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{Read, Write};

/// Header for .oga files.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneHeader {
    pub schema_version: u16,
    pub json_key: String,
    pub name: String,
    pub version: String,
    pub assembly: String,
}

/// In-memory gene annotation database loaded from an .oga file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneIndex {
    pub header: GeneHeader,
    /// Gene symbol -> list of JSON annotations.
    pub genes: HashMap<String, Vec<String>>,
}

impl GeneIndex {
    /// Create a new empty gene index.
    pub fn new(header: GeneHeader) -> Self {
        Self {
            header,
            genes: HashMap::new(),
        }
    }

    /// Add a gene record.
    pub fn add(&mut self, record: GeneRecord) {
        self.genes
            .entry(record.gene_symbol)
            .or_default()
            .push(record.json);
    }

    /// Look up annotations for a gene symbol.
    pub fn get(&self, gene_symbol: &str) -> Option<&Vec<String>> {
        self.genes.get(gene_symbol)
    }

    /// Number of genes in the index.
    pub fn gene_count(&self) -> usize {
        self.genes.len()
    }

    /// Serialize to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(OGA_MAGIC)?;
        writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        let data = bincode::serialize(self)?;
        writer.write_all(&(data.len() as u64).to_le_bytes())?;
        writer.write_all(&data)?;
        Ok(())
    }

    /// Deserialize from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != OGA_MAGIC {
            anyhow::bail!("Invalid OGA magic");
        }
        let mut ver = [0u8; 2];
        reader.read_exact(&mut ver)?;
        if u16::from_le_bytes(ver) != SCHEMA_VERSION {
            anyhow::bail!("Unsupported OGA schema version");
        }
        let mut len_bytes = [0u8; 8];
        reader.read_exact(&mut len_bytes)?;
        let len_u64 = u64::from_le_bytes(len_bytes);
        if len_u64 > MAX_INDEX_PAYLOAD as u64 {
            anyhow::bail!(
                "OGA payload size {} exceeds limit {}",
                len_u64,
                MAX_INDEX_PAYLOAD
            );
        }
        let mut data = vec![0u8; len_u64 as usize];
        reader.read_exact(&mut data)?;
        let index: GeneIndex = bincode::deserialize(&data)?;
        Ok(index)
    }
}

impl GeneAnnotationProvider for GeneIndex {
    fn name(&self) -> &str {
        &self.header.name
    }

    fn json_key(&self) -> &str {
        &self.header.json_key
    }

    fn annotate_gene(&self, gene_symbol: &str) -> Result<Option<String>> {
        match self.genes.get(gene_symbol) {
            Some(jsons) if jsons.len() == 1 => Ok(Some(jsons[0].clone())),
            Some(jsons) => {
                // Multiple annotations: return as JSON array
                let array = format!("[{}]", jsons.join(","));
                Ok(Some(array))
            }
            None => Ok(None),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gene_round_trip() {
        let header = GeneHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "omim".into(),
            name: "OMIM".into(),
            version: "2024-01".into(),
            assembly: "GRCh38".into(),
        };

        let mut index = GeneIndex::new(header);
        index.add(GeneRecord {
            gene_symbol: "BRCA1".into(),
            json: r#"{"mim":113705}"#.into(),
        });
        index.add(GeneRecord {
            gene_symbol: "TP53".into(),
            json: r#"{"mim":191170}"#.into(),
        });

        assert_eq!(index.gene_count(), 2);

        // Serialize and deserialize
        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();
        let loaded = GeneIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();

        assert_eq!(loaded.header.json_key, "omim");
        assert_eq!(loaded.gene_count(), 2);
        assert_eq!(loaded.get("BRCA1").unwrap().len(), 1);
    }

    #[test]
    fn test_gene_provider() {
        let header = GeneHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "omim".into(),
            name: "OMIM".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };

        let mut index = GeneIndex::new(header);
        index.add(GeneRecord {
            gene_symbol: "BRCA1".into(),
            json: r#"{"mim":113705}"#.into(),
        });

        let result = index.annotate_gene("BRCA1").unwrap();
        assert!(result.is_some());
        assert!(result.unwrap().contains("113705"));

        let result = index.annotate_gene("NONEXISTENT").unwrap();
        assert!(result.is_none());
    }
}

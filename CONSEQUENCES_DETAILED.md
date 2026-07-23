# All 49 Consequences with Heuristics

## **RANK 1-7: HIGH IMPACT (Protein severely disrupted)**

### **Rank 1: TranscriptAblation**
**Impact**: HIGH | **Type**: Structural Variant

**What**: Entire transcript completely deleted.

**Heuristic**: Large deletion (>1000bp) that completely covers the transcript span from start to end.

**Example**: 
```
Transcript:    chr1:1000-5000
SV deletion:   chr1:500-6000
               Deletion span > transcript span → TranscriptAblation
```

**Clinical**: Complete loss of gene function. Severe.

---

### **Rank 2: SpliceAcceptorVariant**
**Impact**: HIGH | **Type**: Coding SNV/Indel

**What**: Variant destroys the splice acceptor site (last 2 bases of intron, usually AG).

**Heuristic**: Position is the last 2 bases of an intron (right before exon starts).

**Example**:
```
Exon1[===]Intron[===]AG|Exon2
                    ↑↑ positions 1999-2000 in intron
Variant at 1999 or 2000 → SpliceAcceptorVariant
```

**Clinical**: Splicing broken → no mRNA made. Severe.

---

### **Rank 3: SpliceDonorVariant**
**Impact**: HIGH | **Type**: Coding SNV/Indel

**What**: Variant destroys the splice donor site (first 2 bases of intron, usually GT).

**Heuristic**: Position is the first 2 bases of an intron (right after exon ends).

**Example**:
```
Exon1[===]|GT[===]Intron
          ↑↑ positions 1201-1202 in intron
Variant at 1201 or 1202 → SpliceDonorVariant
```

**Clinical**: Splicing broken → no mRNA made. Severe.

---

### **Rank 4: StopGained**
**Impact**: HIGH | **Type**: Coding SNV/MNV

**What**: SNV/MNV creates a new STOP codon prematurely (mid-sequence, not at end).

**Heuristic**: 
- Original codon: any amino acid
- Mutated codon: TAA, TAG, or TGA (STOP)
- Position: before the normal stop codon

**Example**:
```
Original:  ...Pro-Ala-Gly-Leu-Stop
           ...CCC-GCA-GGA-CTG-TAA
Variant:   ...Pro-Ala-STOP-Leu-Stop
           ...CCC-GCA-TAA-CTG-TAA
                       ↑ new STOP
```

**Clinical**: Truncated protein. Often non-functional. Severe.

---

### **Rank 5: FrameshiftVariant**
**Impact**: HIGH | **Type**: Coding Indel

**What**: Insertion or deletion that shifts the reading frame.

**Heuristic**: 
- Indel length is NOT divisible by 3
- (indel_length % 3) ≠ 0

**Example**:
```
Original:  ATG | GCT | CAA | TAA
           Met | Ala | Gln | STOP

Delete 1 base:  ATG | CTC | AAT | AA? 
           Met | Leu | Asn | ???
               ↑ Frame shifted! All downstream codons wrong
```

**Clinical**: Every codon downstream is corrupted. Severe.

---

### **Rank 6: StopLost**
**Impact**: HIGH | **Type**: Coding SNV/MNV

**What**: The normal STOP codon is changed to an amino acid codon.

**Heuristic**:
- Original codon: TAA, TAG, or TGA (STOP)
- Mutated codon: any amino acid
- Variant causes read-through past natural stop

**Example**:
```
Original:  ...Pro-Ala-Stop
           ...CCC-GCA-TAA
Variant:   ...Pro-Ala-Tyr
           ...CCC-GCA-TAC
                       ↑ changed STOP to Tyr
```

**Clinical**: Protein continues past normal end. Unpredictable function. Severe.

---

### **Rank 7: StartLost**
**Impact**: HIGH | **Type**: Coding SNV/MNV

**What**: The START codon (usually ATG, position 1) is changed.

**Heuristic**:
- Original codon at position 1: ATG (methionine)
- Mutated codon: NOT a start codon (GTG, TTG might be alternate starts, but ATG is canonical)
- Position: codon 1 of CDS

**Example**:
```
Original:  ATG | GCT | CAA | TAA
           Met | Ala | Gln | STOP

Variant:   GTG | GCT | CAA | TAA
           Val | Ala | Gln | STOP
           ↑ not a canonical start
```

**Clinical**: Translation might not initiate. No protein (or wrong start). Severe.

---

### **Rank 8: TranscriptAmplification**
**Impact**: HIGH | **Type**: Structural Variant

**What**: Entire transcript is duplicated.

**Heuristic**: Large duplication (>1000bp) that completely covers the transcript span.

**Example**:
```
Transcript:    chr1:1000-5000
SV duplication: chr1:900-6000
                Duplication span > transcript span → TranscriptAmplification
```

**Clinical**: Gene dosage imbalance. Extra copies of protein. Severe (but less severe than deletion).

---

## **RANK 9-10: HIGH IMPACT (Structural variants)**

### **Rank 9: FeatureElongation**
**Impact**: MODIFIER | **Type**: Structural Variant

**What**: Partial duplication that extends (elongates) the transcript.

**Heuristic**: Duplication partially overlaps transcript but doesn't completely cover it.

**Example**:
```
Transcript:        chr1:1000-5000
SV duplication:    chr1:3000-6000
                   Partial overlap, extends beyond → FeatureElongation
```

**Clinical**: Extra sequence added. Unpredictable effect.

---

### **Rank 10: FeatureTruncation**
**Impact**: MODIFIER | **Type**: Structural Variant

**What**: Partial deletion that truncates (shortens) the transcript.

**Heuristic**: Deletion partially overlaps transcript but doesn't completely cover it.

**Example**:
```
Transcript:   chr1:1000-5000
SV deletion:  chr1:3000-4000
              Partial overlap, shortens → FeatureTruncation
```

**Clinical**: Part of gene deleted. Function likely lost.

---

## **RANK 11-14: MODERATE IMPACT (Protein changed but intact)**

### **Rank 11: InframeInsertion**
**Impact**: MODERATE | **Type**: Coding Indel

**What**: Insertion that preserves the reading frame (length divisible by 3).

**Heuristic**:
- Insertion length % 3 == 0
- Can insert 3, 6, 9, etc. bases

**Example**:
```
Original:  ATG | GCT | CAA | TAA
           Met | Ala | Gln | STOP

Insert 3bp (CAG):  ATG | CAG | GCT | CAA | TAA
           Met | Gln | Ala | Gln | STOP
               ↑ inserted codon
```

**Clinical**: Extra amino acid(s) added. Protein longer. Might disrupt function or be tolerated.

---

### **Rank 12: InframeDeletion**
**Impact**: MODERATE | **Type**: Coding Indel

**What**: Deletion that preserves the reading frame (length divisible by 3).

**Heuristic**:
- Deletion length % 3 == 0
- Can delete 3, 6, 9, etc. bases

**Example**:
```
Original:  ATG | GCT | CAA | TAA
           Met | Ala | Gln | STOP

Delete 3bp (CAA):  ATG | GCT | TAA
           Met | Ala | STOP
               ↑ Gln deleted
```

**Clinical**: Amino acid(s) removed. Protein shorter. Might disrupt function or be tolerated.

---

### **Rank 13: MissenseVariant**
**Impact**: MODERATE | **Type**: Coding SNV/MNV

**What**: SNV/MNV that changes one amino acid to a different one.

**Heuristic**:
- Original AA ≠ Mutated AA
- Neither is STOP
- Single codon affected

**Example**:
```
Original:  GCA (Ala)
Variant:   TCA (Ser)
           Different amino acids → MissenseVariant
```

**Clinical**: Protein sequence changed. Might disrupt function or be tolerated. Most common disease variants.

---

### **Rank 14: ProteinAlteringVariant**
**Impact**: MODERATE | **Type**: Coding (unusual cases)

**What**: A variant that alters the protein in an unusual way (not standard missense/frameshift/stop).

**Heuristic**: Rare. Might be a complex indel or unusual change that affects protein without a standard consequence.

**Example**: A very rare combinatorial change.

**Clinical**: Protein altered. Effect unpredictable.

---

## **RANK 15-22: LOW IMPACT (Minimal or no amino acid change)**

### **Rank 15: SpliceRegionVariant**
**Impact**: LOW | **Type**: Coding SNV/Indel

**What**: Variant near a splice site but not at the critical 2bp sites.

**Heuristic**: 
- Position is 1-3 bases into exon from boundary, OR
- Position is 3-8 bases into intron from boundary
- NOT the first 2bp of intron (SpliceDonorVariant)
- NOT the last 2bp of intron (SpliceAcceptorVariant)

**Example**:
```
Exon[===]|GT???GG[===]Intron
         12345678↑ position 5 or 6 of intron
           → SpliceRegionVariant (not SpliceDonorVariant)
```

**Clinical**: Splicing might be mildly affected. Usually tolerated.

---

### **Rank 16: SpliceDonorFifthBaseVariant**
**Impact**: LOW | **Type**: Coding SNV/Indel

**What**: Variant at position 5 of the donor site (intron).

**Heuristic**: Position is exactly 5 bases into the intron from exon boundary.

**Example**:
```
Exon|GTAAA-Intron
    12345↑ position 5
```

**Clinical**: Weak splice signal. Position 5 is less critical. Usually tolerated.

---

### **Rank 17: SpliceDonorRegionVariant**
**Impact**: LOW | **Type**: Coding SNV/Indel

**What**: Variant in the donor region (positions 3-6 of intron).

**Heuristic**: Position is 3, 4, 5, or 6 bases into the intron from exon.

**Example**:
```
Exon|GTAAGG-Intron
    123456↑ positions 3-6
```

**Clinical**: Donor sequence context. Weakly affects splicing. Usually tolerated.

---

### **Rank 18: SplicePolypyrimidineTractVariant**
**Impact**: LOW | **Type**: Coding SNV/Indel

**What**: Variant in the polypyrimidine tract (3-17 bases upstream of acceptor site).

**Heuristic**: Position is 3-17 bases before the end of the intron (before acceptor AG).

**Example**:
```
Intron[===]CCCCCCCCC|AG Exon
      ← 17bp to 3bp ↑
      polypyrimidine tract
```

**Clinical**: Helps position spliceosome. Weakly affects splicing. Usually tolerated.

---

### **Rank 19: IncompleteTerminalCodonVariant**
**Impact**: LOW | **Type**: Coding SNV/Indel

**What**: Variant in the terminal (stop) codon region, but not a complete stop loss.

**Heuristic**: Rare. Codon at end of CDS is partially affected but still acts as stop.

**Example**: Last base of stop codon mutated but stop still functions.

**Clinical**: Usually tolerated.

---

### **Rank 20: StartRetainedVariant**
**Impact**: LOW | **Type**: Coding SNV

**What**: Variant in the start codon but a different start codon is still functional.

**Heuristic**: 
- Original codon: ATG (canonical start)
- Mutated codon: GTG or TTG (alternate start codons that are still valid)
- Position: codon 1

**Example**:
```
Original: ATG (Met, standard start)
Variant:  GTG (Val, but also a valid start codon in some genes)
```

**Clinical**: Translation still initiates (just different codon). Usually tolerated.

---

### **Rank 21: StopRetainedVariant**
**Impact**: LOW | **Type**: Coding SNV

**What**: Variant in the stop codon but it's still a stop codon (different stop codon).

**Heuristic**:
- Original codon: TAA, TAG, or TGA (any STOP)
- Mutated codon: TAA, TAG, or TGA (still a STOP, just different)

**Example**:
```
Original: TAA (stop)
Variant:  TAG (also stop)
```

**Clinical**: Still stops translation. Usually harmless.

---

### **Rank 22: SynonymousVariant**
**Impact**: LOW | **Type**: Coding SNV/MNV

**What**: SNV/MNV that doesn't change the amino acid (silent mutation).

**Heuristic**:
- Original codon AA == Mutated codon AA
- "Wobble" base (3rd position) often involved

**Example**:
```
Original: GCA (Ala)
Variant:  GCT (Ala)
          Same amino acid → SynonymousVariant
```

**Clinical**: Protein unchanged. Usually harmless (though can affect splicing/codon usage rarely).

---

## **RANK 23-31: MODIFIER IMPACT (Non-coding or catch-all)**

### **Rank 23: CodingSequenceVariant**
**Impact**: MODIFIER | **Type**: Coding (catch-all)

**What**: Variant in a coding sequence, but consequence couldn't be determined precisely.

**Heuristic**: Variant is in CDS but doesn't match any specific consequence (fallback).

**Example**: Edge case where prediction failed.

**Clinical**: In coding region but effect unknown.

---

### **Rank 24: MatureMirnaVariant**
**Impact**: MODIFIER | **Type**: Non-coding RNA

**What**: Variant within a mature miRNA sequence.

**Heuristic**: Variant position is within the miRNA exon (20-25bp sequence).

**Example**:
```
miRNA: 5' UAGCUUAUCAGACUGAUGUUGA 3'
              ↑ variant here
```

**Clinical**: Affects target binding. Effect depends on position within miRNA.

---

### **Rank 25: FivePrimeUtrVariant**
**Impact**: MODIFIER | **Type**: Coding RNA (untranslated)

**What**: Variant in the 5' UTR (before start codon).

**Heuristic**: 
- In exon: YES
- Before CDS start: YES (forward strand) OR after CDS end (reverse strand)
- In coding transcript: YES

**Example**:
```
[5'UTR][CDS][3'UTR]
 ↑
Variant here
```

**Clinical**: Transcribed but not translated. Usually harmless (unless affects Kozak sequence near start).

---

### **Rank 26: ThreePrimeUtrVariant**
**Impact**: MODIFIER | **Type**: Coding RNA (untranslated)

**What**: Variant in the 3' UTR (after stop codon).

**Heuristic**:
- In exon: YES
- After CDS end: YES (forward strand) OR before CDS start (reverse strand)
- In coding transcript: YES

**Example**:
```
[5'UTR][CDS][3'UTR]
              ↑
           Variant here
```

**Clinical**: Transcribed but not translated. Usually harmless (unless affects polyadenylation or mRNA stability).

---

### **Rank 27: NonCodingTranscriptExonVariant**
**Impact**: MODIFIER | **Type**: Non-coding transcript

**What**: Variant in an exon of a non-coding transcript (lncRNA, miRNA, etc.).

**Heuristic**:
- Transcript biotype: NOT protein_coding
- In exon: YES

**Example**: Variant in lncRNA exon.

**Clinical**: Affects non-coding RNA. Usually harmless.

---

### **Rank 28: IntronVariant**
**Impact**: MODIFIER | **Type**: Any transcript

**What**: Variant within an intron, away from splice sites.

**Heuristic**:
- In intron: YES
- Not at splice donor/acceptor: YES
- >8 bases from exon boundaries: YES

**Example**:
```
Exon[===]IntronVARIANT_HERE[===]Exon
            ↑ 50bp into intron
```

**Clinical**: Intron removed during splicing. Usually harmless (unless affects splicing signals).

---

### **Rank 29: NmdTranscriptVariant**
**Impact**: MODIFIER | **Type**: Modifier flag

**What**: Transcript is marked as nonsense_mediated_decay biotype (unstable).

**Heuristic**: Check `transcript.biotype == "nonsense_mediated_decay"`

**Example**: BRCA1-002 isoform.

**Clinical**: Already unstable transcript. Variant here compounds the instability.

---

### **Rank 30: NonCodingTranscriptVariant**
**Impact**: MODIFIER | **Type**: Non-coding transcript

**What**: Generic consequence for non-coding transcript (fallback).

**Heuristic**: 
- Transcript biotype: lncRNA, miRNA, etc.
- Consequence couldn't be more specific

**Example**: Variant in ncRNA without being in a specific exon.

**Clinical**: Non-coding variant. Usually harmless.

---

### **Rank 31: CodingTranscriptVariant**
**Impact**: MODIFIER | **Type**: Coding transcript (catch-all)

**What**: Variant in coding transcript but consequence unclear (fallback).

**Heuristic**: Coding transcript but specific consequence couldn't be determined.

**Example**: Edge case or unusual situation.

**Clinical**: In coding region but effect unclear.

---

## **RANK 32-33: UPSTREAM/DOWNSTREAM**

### **Rank 32: UpstreamGeneVariant**
**Impact**: MODIFIER | **Type**: Distance-based

**What**: Variant is upstream of the gene, within distance threshold (default 5000bp).

**Heuristic**:
- Does NOT overlap transcript: YES
- Distance from transcript: < 5000bp
- Direction: UPSTREAM (before gene start on forward strand, after on reverse)

**Example**:
```
Variant──────2000bp──────[Gene Start]
         upstream
```

**Clinical**: Might affect promoter or regulatory region. Usually harmless.

---

### **Rank 33: DownstreamGeneVariant**
**Impact**: MODIFIER | **Type**: Distance-based

**What**: Variant is downstream of the gene, within distance threshold (default 5000bp).

**Heuristic**:
- Does NOT overlap transcript: YES
- Distance from transcript: < 5000bp
- Direction: DOWNSTREAM (after gene end on forward strand, before on reverse)

**Example**:
```
[Gene End]──────2000bp──────Variant
        downstream
```

**Clinical**: Might affect regulatory region or polyadenylation. Usually harmless.

---

## **RANK 34-39: REGULATORY REGION**

### **Rank 34: TfbsAblation**
**Impact**: HIGH | **Type**: Structural Variant

**What**: Deletion completely removes a transcription factor binding site.

**Heuristic**: SV deletion completely covers a TFBS region.

**Example**: Deletion of known TF binding site.

**Clinical**: Loss of transcription regulation. Severe for that TF.

---

### **Rank 35: TfbsAmplification**
**Impact**: MODERATE | **Type**: Structural Variant

**What**: Duplication of a transcription factor binding site.

**Heuristic**: SV duplication affects a TFBS region.

**Example**: Extra copies of TF binding site.

**Clinical**: Altered transcription regulation. Moderate effect.

---

### **Rank 36: TfBindingSiteVariant**
**Impact**: MODIFIER | **Type**: SNV/Indel

**What**: Variant within a known transcription factor binding site.

**Heuristic**: Position overlaps annotated TFBS in genome database.

**Example**: SNV in a known TF binding site sequence.

**Clinical**: Might affect TF binding. Usually weak effect.

---

### **Rank 37: RegulatoryRegionAblation**
**Impact**: HIGH | **Type**: Structural Variant

**What**: Deletion completely removes a regulatory region.

**Heuristic**: SV deletion completely covers an annotated regulatory region (promoter, enhancer, etc.).

**Example**: Large deletion removing promoter region.

**Clinical**: Loss of gene regulation. Severe.

---

### **Rank 38: RegulatoryRegionAmplification**
**Impact**: MODERATE | **Type**: Structural Variant

**What**: Duplication of a regulatory region.

**Heuristic**: SV duplication affects annotated regulatory region.

**Example**: Extra copies of enhancer.

**Clinical**: Altered regulation. Moderate effect.

---

### **Rank 39: RegulatoryRegionVariant**
**Impact**: MODIFIER | **Type**: SNV/Indel

**What**: Variant within a known regulatory region (promoter, enhancer, silencer).

**Heuristic**: Position overlaps annotated regulatory region.

**Example**: SNV in promoter region.

**Clinical**: Might affect regulation. Usually weak effect.

---

## **RANK 40-41: INTERGENIC / GENERIC**

### **Rank 40: IntergenicVariant**
**Impact**: MODIFIER | **Type**: Distance-based

**What**: Variant is between genes, far from any gene (>5000bp away).

**Heuristic**:
- Does NOT overlap any transcript: YES
- Distance to nearest gene: > 5000bp

**Example**:
```
[Gene A]────────50000bp────────[Gene B]
                  ↑
              IntergenicVariant
```

**Clinical**: Between genes. Usually harmless.

---

### **Rank 41: SequenceVariant**
**Impact**: MODIFIER | **Type**: Generic catch-all

**What**: A variant detected, but no specific consequence could be determined.

**Heuristic**: Fallback when nothing else applies.

**Example**: Unknown variant type.

**Clinical**: Unknown effect.

---

## **RANK 42-47: STRUCTURAL VARIANTS - Copy Number**

### **Rank 42: CopyNumberChange**
**Impact**: MODIFIER | **Type**: SV

**What**: Copy number varies (could be gain or loss, indeterminate).

**Heuristic**: CNV without clear direction (e.g., copy number polymorphism).

**Example**: Region with variable copy number.

**Clinical**: Effect depends on size and direction.

---

### **Rank 43: CopyNumberIncrease**
**Impact**: MODIFIER | **Type**: SV

**What**: Duplication or copy number gain (more copies than normal).

**Heuristic**: SV type = Duplication or CopyNumberGain AND partial or full overlap.

**Example**: 3 copies of region instead of 2.

**Clinical**: Gene dosage imbalance (too much protein). Often pathogenic.

---

### **Rank 44: CopyNumberDecrease**
**Impact**: MODIFIER | **Type**: SV

**What**: Deletion or copy number loss (fewer copies than normal).

**Heuristic**: SV type = Deletion or CopyNumberLoss AND partial or full overlap.

**Example**: 1 copy of region instead of 2 (monosomy).

**Clinical**: Gene dosage imbalance (too little protein). Often pathogenic.

---

## **RANK 45-47: STRUCTURAL VARIANTS - STR**

### **Rank 45: ShortTandemRepeatChange**
**Impact**: MODIFIER | **Type**: SV

**What**: Copy number of a short tandem repeat (STR) changes indeterminately.

**Heuristic**: STR region with changed repeat count, direction unclear.

**Example**: Microsatellite with variable repeats.

**Clinical**: Variable effect based on repeat count.

---

### **Rank 46: ShortTandemRepeatExpansion**
**Impact**: MODIFIER | **Type**: SV

**What**: Short tandem repeat expands (more repeats).

**Heuristic**: STR region with increased repeat count.

**Example**: CAG repeats in Huntingtin gene (normal: 15-35, disease: >40).

**Clinical**: Often pathogenic. Anticipation (worse in next generation).

---

### **Rank 47: ShortTandemRepeatContraction**
**Impact**: MODIFIER | **Type**: SV

**What**: Short tandem repeat contracts (fewer repeats).

**Heuristic**: STR region with decreased repeat count.

**Example**: GAA repeats in Frataxin gene contracting.

**Clinical**: Might restore function if originally expanded.

---

## **RANK 48: GENE FUSIONS**

### **Rank 48: UnidirectionalGeneFusion**
**Impact**: MODIFIER | **Type**: SV

**What**: Two genes fused together, one gene's sequence inserted into another.

**Heuristic**: SV breakend (translocation) that joins two different genes.

**Example**: BCR-ABL fusion (Philadelphia chromosome in chronic myeloid leukemia).

**Clinical**: Fusion protein with new function. Often oncogenic or loss of function for one gene.

---

## **RANK 49: GENERIC SV**

### **Rank 49: TranscriptVariant**
**Impact**: MODIFIER | **Type**: SV catch-all

**What**: Structural variant overlaps transcript, but specific consequence indeterminate.

**Heuristic**: SV overlaps transcript but type/impact unclear.

**Example**: Inversion or complex SV affecting transcript.

**Clinical**: Transcript affected but effect unclear.

---

## **Quick Reference by Common Scenarios**

| Scenario | Consequence |
|----------|-------------|
| SNV changes amino acid | MissenseVariant |
| SNV doesn't change amino acid | SynonymousVariant |
| SNV creates STOP mid-sequence | StopGained |
| 1bp deletion in CDS | FrameshiftVariant |
| 3bp deletion in CDS | InframeDeletion |
| Variant at first 2bp of intron | SpliceDonorVariant |
| Variant at last 2bp of intron | SpliceAcceptorVariant |
| Variant 50bp into intron | IntronVariant |
| Variant in 5' UTR | FivePrimeUtrVariant |
| Variant in 3' UTR | ThreePrimeUtrVariant |
| Variant 2000bp before gene | UpstreamGeneVariant |
| Variant 2000bp after gene | DownstreamGeneVariant |
| Variant 50000bp away | IntergenicVariant |
| Large deletion covering transcript | TranscriptAblation |
| Large duplication covering transcript | TranscriptAmplification |
| Variant in lncRNA exon | NonCodingTranscriptExonVariant |


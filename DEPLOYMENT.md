# Deploying fastVEP on Hetzner Cloud

Complete guide to hosting fastVEP as a public web service with full human genome annotation (GRCh38) and all supplementary databases.

[fastvep.org](https://fastvep.org) runs this guide on one server with 4 vCPU, 16 GB of RAM and 160 GB of disk.
It serves nine genomes from one process, with nginx in front.

## Cost Overview

| Component | Monthly Cost |
|-----------|-------------|
| Hetzner CCX23 (4 vCPU, 16GB RAM, 160GB SSD) | ~$30 |
| Cloudflare DNS + CDN (free tier) | $0 |
| **Total** | **~$30/month** |

Upgrade to CCX33 (8 vCPU, 32GB RAM, ~$55/month) if you want headroom for heavy traffic.

## Deployment Tiers

fastVEP ships with ready-to-run deployment scripts in `scripts/`. Pick the one that fits your server and use case.

| Tier | Script | Disk used | What you get |
|:-----|:-------|:----------|:-------------|
| **Minimal** | `scripts/deploy-minimal.sh` | ~5 GB | Human GRCh38 + ClinVar |
| **Clinical** | `scripts/deploy-clinical.sh` | ~23 GB | GRCh38 + ClinVar + dbSNP + seven other organisms (mouse, rat, zebrafish, fruit fly, C. elegans, yeast, Arabidopsis) |
| **Full Research** | `scripts/deploy-full.sh` | ~26 GB + gnomAD | GRCh38 + GRCh37 + ClinVar (both) + dbSNP (GRCh38) + gnomAD v4.1 genomes (GRCh38) + the same seven organisms |

The disk figures are what the files take on fastvep.org (8 November 2026): the whole data directory, nine genomes with dbSNP and ClinVar on GRCh38 and ClinVar on GRCh37, is 26 GB, of which GRCh38 is 14 GB and its dbSNP 9.4 GB.
gnomAD is not in that figure: see the sizes below.

### Quick deploy with scripts

After building the binaries (Section 2), run the deployment script **on the server**:

```bash
# 1. SSH into the server
ssh -i ~/.ssh/fastvep_key root@Your.server.IP

# 2. Make sure the repo is up-to-date (scripts are in the repo)
cd /root/fastVEP
git pull

# 3. Rebuild both binaries (CLI + web server)
cargo build --release

# 4. Copy binaries to /opt/fastvep/bin
mkdir -p /opt/fastvep/bin
cp target/release/fastvep-web /opt/fastvep/bin/fastvep-web
cp target/release/fastvep     /opt/fastvep/bin/fastvep
chmod +x /opt/fastvep/bin/fastvep-web /opt/fastvep/bin/fastvep

# 5. Run the deployment script (pick one)
export FASTVEP_BIN=/root/fastVEP/target/release
bash scripts/deploy-full.sh        # everything the site serves
# OR: bash scripts/deploy-clinical.sh   # no GRCh37, no gnomAD
# OR: bash scripts/deploy-minimal.sh    # quickest, just hg38 + ClinVar
```

**Alternatively**, if you built locally (Mac) and want to upload the scripts to the server:

```bash
# From your local machine
scp scripts/deploy-full.sh root@Your.server.IP:/root/fastVEP/scripts/
ssh -i ~/.ssh/fastvep_key root@Your.server.IP 'cd /root/fastVEP && bash scripts/deploy-full.sh'
```

The scripts are **idempotent**: you can re-run them safely if interrupted.
They skip a file or a database that already exists, and gnomAD is built one chromosome at a time, each into a database of its own, so an interrupted build resumes at the chromosome it stopped in.

> **Which tier should I pick?**
> - **Minimal** (~5 GB): fast variant consequence + ClinVar.
> - **Clinical** (~23 GB once built): adds dbSNP rs IDs and the seven other organisms. The dbSNP build needs 29.6 GB of download space while it runs, and the script deletes the download afterwards.
> - **Full Research** (~26 GB plus gnomAD, several hours): adds GRCh37 and gnomAD v4.1 population frequencies.
>
> All three fit a 160 GB disk.
> The scripts do not build the rest of the ACMG reference stack: see [ACMG on the web needs the reference stack](#acmg-on-the-web-needs-the-reference-stack).

You can always start with Minimal and upgrade later: the scripts build on each other.

> **Disk space warning**: the gnomAD v4.1 genomes VCFs are 0.6 to 47 GB **per chromosome** (563 GB in all, chr2 is the largest) and the exomes 0.1 to 19 GB (198 GB in all, chr1 is the largest; measured on 8 November 2026).
> The Full Research script downloads one chromosome at a time, builds its database, then deletes the VCF, so the temporary space it needs is the size of the largest VCF: 47 GB for the genomes.

## 1. Provision the Server

### Create account and server

1. Sign up at [hetzner.com/cloud](https://www.hetzner.com/cloud)
2. Create a new project
3. Add an SSH key (Settings → SSH Keys)
4. Create a server:
   - **Location**: Ashburn (us-east) or Falkenstein (eu-central) — pick closest to your users
   - **Image**: Ubuntu 24.04
   - **Type**: CCX23 (Dedicated, 4 vCPU AMD, 16GB RAM, 80GB SSD)
   - **SSH key**: select yours
   - **Name**: `fastvep`

### Domain Setup (Hostinger)

To connect your domain from Hostinger to this server:

1. Log in to your **Hostinger Control Panel**.
2. Go to **Domains** → Select your domain → **DNS / Nameservers**.
3. Add or edit an **A Record**:
   - **Type**: `A`
   - **Name**: `@` (for the root domain) or `fastvep` (for a subdomain like `fastvep.yourdomain.com`).
   - **Points to**: `Your.server.IP` (Your server IP).
   - **TTL**: Default (usually 14400 or 3600).
4. Save the record. Propagation may take up to 24 hours but is usually faster.

### Initial setup

```bash
ssh -i ~/.ssh/fastvep_key root@Your.server.IP

# System updates
apt update && apt upgrade -y

# Create service user
useradd -m -s /bin/bash fastvep
mkdir -p /opt/fastvep/data /opt/fastvep/bin
chown -R fastvep:fastvep /opt/fastvep

# Install nginx
apt install -y nginx certbot python3-certbot-nginx

# Firewall
ufw allow OpenSSH
ufw allow 'Nginx Full'
ufw enable
```

Then switch SSH to keys only.
sshd takes the first value it reads for a keyword, and it reads the files of `/etc/ssh/sshd_config.d/` in name order before the rest of `sshd_config`.
A file named `99-...` therefore cannot undo a `PasswordAuthentication yes` from an earlier one (checked on Ubuntu 24.04 with a competing `50-` file), so name yours to be read first, and look at the effective values instead of trusting the file:

```bash
printf 'PasswordAuthentication no\nPermitRootLogin prohibit-password\nMaxAuthTries 3\n' > /etc/ssh/sshd_config.d/00-fastvep.conf
sshd -t && systemctl reload ssh
sshd -T | grep -E '^(passwordauthentication|permitrootlogin|maxauthtries) '
```

Keep the session you are in open, and log in again from a second terminal with your key before you close it.

## 2. Prepare and Build the Binary

On your development machine, you'll need the source code and Rust installed.

### Which build method should I use?

| Method | Pros | Cons |
|:---|:---|:---|
| **Build on Server** | **Simplest.** No local setup. Guaranteed compatibility. | Uses server CPU/RAM (may slow down site during build). |
| **Cross-compile** | Faster (uses your Mac/PC). No extra tools on server. | Requires Docker/Colima locally. More complex to set up. |

**Recommendation**: If it's your first time, **Build on the server directly**. It's much harder to mess up.

### Clone the repository
```bash
git clone https://github.com/Huang-lab/fastVEP.git
cd fastVEP
```

### Build the binary
If building on a Mac/Windows for a Linux server (cross-compiling):

```bash
# 1. Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# 2. Install build dependencies (required for compiling Rust tools)
# On Ubuntu/Debian:
sudo apt update && sudo apt install -y build-essential pkg-config libssl-dev
# On Mac:
# xcode-select --install

# 3. Install cross for cross-compilation
cargo install cross --git https://github.com/cross-rs/cross

# 4. Build for Linux (requires Docker or Colima)
# Note: If building on the server, you must install Docker:
# sudo apt install -y docker.io && sudo usermod -aG docker $USER && newgrp docker
cross build --release --target x86_64-unknown-linux-gnu -p fastvep-web

# Upload binary
scp target/x86_64-unknown-linux-gnu/release/fastvep-web root@YOUR_SERVER_IP:/opt/fastvep/bin/
ssh root@YOUR_SERVER_IP 'chmod +x /opt/fastvep/bin/fastvep-web'
```

**Alternative: Build on the server directly**
This is often easier if you don't want to set up Docker/cross locally.
```bash
ssh root@YOUR_SERVER_IP
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Install build dependencies
apt update && apt install -y build-essential pkg-config libssl-dev

# Clone and build
git clone https://github.com/Huang-lab/fastVEP.git
cd fastVEP
cargo build --release
mkdir -p /opt/fastvep/bin
cp target/release/fastvep-web /opt/fastvep/bin/fastvep-web
cp target/release/fastvep     /opt/fastvep/bin/fastvep
chmod +x /opt/fastvep/bin/fastvep-web /opt/fastvep/bin/fastvep
```

> **Note**: `cargo build --release` (without `-p`) builds both the `fastvep-web` server and the `fastvep` CLI. The CLI is needed to build SA databases (ClinVar, gnomAD, dbSNP, etc.).

## 3. Download Genome Data

fastVEP supports multiple genomes via the `--data-dir` flag. Each subdirectory is an independent genome (different species or different human reference builds). The web UI shows them in a dropdown labeled with the build/source so users pick unambiguously.

### Directory layout

```
/opt/fastvep/data/
├── hg38_ensembl_115/              # ← each subdir = one genome
│   ├── Homo_sapiens.GRCh38.115.gff3
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
│   └── sa/
│       ├── clinvar.osa2
│       ├── gnomad_chr1.osa2     # one file per chromosome
│       └── ...
├── hg37_ensembl_87/
│   ├── Homo_sapiens.GRCh37.87.gff3
│   └── ...
├── mouse_grcm39_ensembl_115/
│   ├── Mus_musculus.GRCm39.115.gff3
│   └── ...
├── arabidopsis_tair10_ensemblplants_63/
│   ├── Arabidopsis_thaliana.TAIR10.63.gff3
│   └── ...
└── ...
```

The web page builds the dropdown from the directory names, so **name them for what is in them**: the species or build, the annotation source and its release (`hg38_ensembl_115`, `hg37_ensembl_87`, `arabidopsis_tair10_ensemblplants_63`), never `genome1`.
A name that says a release the files are not from is worse than a vague one: fastvep.org once had a GRCh37 directory named for release 75 that held the release 87 gene set.
The web server serves every subdirectory that holds a GFF3, and a request names the genome by its directory name.

SSH into the server:

```bash
ssh root@YOUR_SERVER_IP
su - fastvep
mkdir -p /opt/fastvep/data
cd /opt/fastvep/data
```

If you built on the server (recommended), both `fastvep` and `fastvep-web` are already in `/opt/fastvep/bin/`. If you cross-compiled, also upload the CLI binary:

```bash
# Only needed if you cross-compiled — skip if you built on the server
scp target/x86_64-unknown-linux-gnu/release/fastvep root@YOUR_SERVER_IP:/opt/fastvep/bin/
chmod +x /opt/fastvep/bin/fastvep
```

Install `samtools` for FASTA indexing:

```bash
apt install -y samtools
```

---

### Human — GRCh38 / hg38 (Ensembl 115, primary assembly)

**This is the recommended default.** Most modern VCFs are aligned to GRCh38.

```bash
mkdir -p hg38_ensembl_115/sa && cd hg38_ensembl_115

# Gene annotation (~1.2GB uncompressed)
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
gunzip Homo_sapiens.GRCh38.115.gff3.gz

# Reference FASTA (~850MB compressed → ~2.9GB uncompressed)
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

cd /opt/fastvep/data
```

### Human - GRCh37 / hg37 (Ensembl 87 gene set)

Still widely used for clinical and legacy data.
The Ensembl gene set of GRCh37 is release 87.
`Homo_sapiens.GRCh37.87.gff3.gz` is the same file (37,903,989 bytes, 20 March 2017) under `release-87`, `release-113` and `current` of `ftp.ensembl.org/pub/grch37/`, so the directory is named after the gene set and not after a release number that changes nothing in it.

```bash
mkdir -p hg37_ensembl_87/sa && cd hg37_ensembl_87

wget https://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
gunzip Homo_sapiens.GRCh37.87.gff3.gz

wget https://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa

cd /opt/fastvep/data
```

### Human — CHM13 / T2T (GCA_009914755.4)

Telomere-to-telomere complete reference. Use if your VCFs are aligned to T2T-CHM13v2.0.

```bash
mkdir -p t2t_chm13v2/sa && cd t2t_chm13v2

# Gene annotation (CAT+Liftoff GFF3)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/vertebrate_mammalian/Homo_sapiens/T2T-CHM13v2.0/GCF_009914755.1-RS_2023_10/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
mv GCF_009914755.1_T2T-CHM13v2.0_genomic.gff t2t_chm13v2.gff3

# Reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
mv GCA_009914755.4_T2T-CHM13v2.0_genomic.fna t2t_chm13v2.fa
samtools faidx t2t_chm13v2.fa

cd /opt/fastvep/data
```

### Human — GRCh38 (GENCODE 47 alternative to Ensembl)

If you prefer GENCODE's annotation over Ensembl's (they're nearly identical but numbered differently):

```bash
mkdir -p hg38_gencode_47/sa && cd hg38_gencode_47

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz
gunzip gencode.v47.annotation.gff3.gz

# FASTA: reuse the Ensembl primary_assembly FASTA (same sequence), or download GENCODE's:
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa

cd /opt/fastvep/data
```

---

### Model organisms

All but Arabidopsis are available from Ensembl release 115.
Same pattern: `gff3/<species>/` and `fasta/<species>/dna/`.

```bash
cd /opt/fastvep/data

# Mouse (Mus musculus) — GRCm39
mkdir -p mouse_grcm39_ensembl_115 && cd mouse_grcm39_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz
gunzip Mus_musculus.GRCm39.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
cd /opt/fastvep/data

# Rat (Rattus norvegicus) - GRCr8
mkdir -p rat_grcr8_ensembl_115 && cd rat_grcr8_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gff3.gz
gunzip Rattus_norvegicus.GRCr8.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz
gunzip Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz
samtools faidx Rattus_norvegicus.GRCr8.dna.toplevel.fa
cd /opt/fastvep/data

# Zebrafish (Danio rerio) — GRCz11
mkdir -p zebrafish_grcz11_ensembl_115 && cd zebrafish_grcz11_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/danio_rerio/Danio_rerio.GRCz11.115.gff3.gz
gunzip Danio_rerio.GRCz11.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
samtools faidx Danio_rerio.GRCz11.dna.primary_assembly.fa
cd /opt/fastvep/data

# Fruit fly (Drosophila melanogaster) - BDGP6.54
mkdir -p drosophila_bdgp6_ensembl_115 && cd drosophila_bdgp6_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gff3.gz
gunzip Drosophila_melanogaster.BDGP6.54.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz
samtools faidx Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa
cd /opt/fastvep/data

# C. elegans — WBcel235
mkdir -p celegans_wbcel235_ensembl_115 && cd celegans_wbcel235_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz
gunzip Caenorhabditis_elegans.WBcel235.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
samtools faidx Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
cd /opt/fastvep/data

# Yeast (S. cerevisiae) — R64-1-1
mkdir -p yeast_r64_ensembl_115 && cd yeast_r64_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
cd /opt/fastvep/data

# Arabidopsis - TAIR10 (Ensembl Plants has its own release numbers: 63 here)
mkdir -p arabidopsis_tair10_ensemblplants_63 && cd arabidopsis_tair10_ensemblplants_63
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3.gz
gunzip Arabidopsis_thaliana.TAIR10.63.gff3.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
cd /opt/fastvep/data
```

Ensembl Plants numbers its releases separately from Ensembl's (`pub/plants/release-63` beside `pub/release-115`), so the directory name says which one it is.
The TAIR10 sequence file is byte-identical in releases 58 and 63 (SHA-256 checked on 8 November 2026), so moving a server from one to the other changes the gene models and not the FASTA.
`ftp.ensemblgenomes.org` answers on plain HTTP only; use the `ebi.ac.uk` host above for HTTPS.

---

### Build supplementary annotation databases (human)

SA databases are build-specific.
ClinVar, gnomAD and dbSNP have separate GRCh37 and GRCh38 files: use the right one for each genome subdirectory, and keep one `sa/` directory to one build.

`sa-build` writes the database under the name it is given, and it replaces a file of that name.
A source that comes as one file per chromosome therefore needs a name per chromosome (`gnomad_chr1`, `gnomad_chr2`, and so on).
One name for all of them leaves the database of the last chromosome alone, which says "no record" for the other twenty-two (checked with dbSNP-format input: building two inputs into one name leaves a file byte-identical to the second input built alone).
A build that finds no record is an error (`--allow-empty` overrides it), because an empty database answers "no record" to every query, and for a frequency source that reads as "absent from the population".

For **hg38_ensembl_115** (GRCh38):

```bash
cd /opt/fastvep/data/hg38_ensembl_115/sa

# ClinVar (about 200 MB download, 43 MB .osa2)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar \
  -i clinvar.vcf.gz -o clinvar --assembly GRCh38
rm clinvar.vcf.gz

# dbSNP (29.6 GB download, 9.4 GB .osa2)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
# --chunk-bits 16: dbSNP is dense (about 3 MB of records in each megabase), and at the default width a scattered
# lookup costs 23 ms, so a VCF of 1,000 variants spread over the genome took 24 s (see docs/SUPPLEMENTARY_ANNOTATIONS.md)
/opt/fastvep/bin/fastvep sa-build --source dbsnp \
  -i GCF_000001405.40.gz -o dbsnp --assembly GRCh38 --chunk-bits 16
rm GCF_000001405.40.gz

# gnomAD v4.1 genomes (0.6-47 GB per chromosome VCF, 563 GB in all)
# IMPORTANT: download ONE chromosome at a time, build its database, then delete the VCF.
for chr in {1..22} X; do
  VCF="gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
  wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/${VCF}"
  /opt/fastvep/bin/fastvep sa-build --source gnomad \
    -i "$VCF" -o "gnomad_chr${chr}" --assembly GRCh38 --chunk-bits 16   # a name per chromosome
  rm -f "$VCF"   # delete raw VCF to free disk space
  echo "chr${chr} done. Disk free: $(df -h --output=avail /opt/fastvep/data | tail -1 | xargs)"
done
```

The gnomAD exomes (`.../release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz`, 0.1-19 GB per chromosome, 198 GB in all) are the smaller release, and the one the ACMG benchmark reads; [docs/ACMG_SETUP.md](docs/ACMG_SETUP.md) says which release to build for classification and why.
A pilot on the fastvep.org server built chromosome 21 of the exomes: 2.25 GB of VCF in 153 s on one core, 2,188,842 records, a 20.7 MB `.osa2`, 89 MB of memory at its peak.
Scaled by download size that is about 193 million records, a 1.8 GB database and 3.7 core-hours for all chromosomes (an extrapolation from one chromosome, not a measurement).
The build used one core (99% CPU), so chromosomes can be built side by side, one per core.

The other sources (1000 Genomes, REVEL, AlphaMissense, PhyloP, SpliceAI, gnomAD gene constraints, ClinGen, RepeatMasker) are built by [`scripts/build-sa-databases.sh`](scripts/build-sa-databases.sh).
[docs/SA_DATABASES.md](docs/SA_DATABASES.md) lists every source, what its download costs, and which need an account.

For **hg37_ensembl_87** (GRCh37):

```bash
cd /opt/fastvep/data/hg37_ensembl_87/sa

# ClinVar GRCh37
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar \
  -i clinvar.vcf.gz -o clinvar --assembly GRCh37
rm clinvar.vcf.gz

# dbSNP GRCh37 (28.2 GB download; dbSNP build 157, GRCh37.p13, contigs named NC_000001.10 and so on).
# sa-build maps those contigs for --assembly GRCh37 and refuses to write a database when none matched.
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
/opt/fastvep/bin/fastvep sa-build --source dbsnp \
  -i GCF_000001405.25.gz -o dbsnp --assembly GRCh37 --chunk-bits 16
rm GCF_000001405.25.gz

# gnomAD v2.1.1 for GRCh37 (v3 and v4 are GRCh38 only)
# Same pattern: one chromosome at a time, a name per chromosome, delete the VCF.
# gnomAD v2.1.1 has no Y chromosome file; that is expected.
for chr in {1..22} X; do
  VCF="gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz"
  wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/${VCF}"
  /opt/fastvep/bin/fastvep sa-build --source gnomad \
    -i "$VCF" -o "gnomad_chr${chr}" --assembly GRCh37 --chunk-bits 16
  rm -f "$VCF"
done
```

---

### Disk space estimate per genome

Human GRCh38, as it is on fastvep.org (8 November 2026):

| File | Size |
|------|------|
| GFF3 | 1.2 GB |
| GFF3 cache (built on the first start) | 292 MB |
| FASTA + .fai | 3.2 GB |
| ClinVar .osa2 | 43 MB |
| dbSNP .osa2 | 9.4 GB |
| gnomAD exomes .osa2, all chromosomes | about 1.8 GB (extrapolated from chromosome 21, 20.7 MB) |
| gnomAD genomes .osa2 | not measured |
| The other ACMG sources | see [docs/ACMG_SETUP.md](docs/ACMG_SETUP.md) |

The whole data directory of fastvep.org, nine genomes, is 26 GB: GRCh38 14 GB, GRCh37 3.6 GB, mouse 3.2 GB, rat 3.1 GB, zebrafish 1.6 GB, and the four small genomes (fruit fly, C. elegans, Arabidopsis, yeast) 0.7 GB together.
The gnomAD build needs, besides its database, the size of the largest VCF as temporary space (47 GB for the genomes, 19 GB for the exomes).

> **Disk planning**: the nine genomes with ClinVar and dbSNP, the gnomAD build's temporary space and the databases of the ACMG stack fit on the 160 GB disk of a CCX23.
> Attach a Hetzner Volume instead if you add genomes or keep the raw VCFs.

## 4. Configure the Service

### systemd unit

One unit serves every genome of the data directory: a request names the genome it wants, nothing is replaced when a visitor chooses, and so `--allow-model-replacement` stays off.
The service runs as the `fastvep` user and not as root, so the data directory has to belong to that user (the server writes its transcript caches beside the GFF3 files), and its counters go to a state directory of their own.

```bash
chown -R fastvep:fastvep /opt/fastvep/data
mkdir -p /opt/fastvep/state && chown fastvep:fastvep /opt/fastvep/state

cat > /etc/systemd/system/fastvep.service << 'EOF'
[Unit]
Description=fastVEP Variant Annotation Server
Documentation=https://github.com/Huang-lab/fastVEP/blob/master/DEPLOYMENT.md
After=network.target

[Service]
Type=simple
User=fastvep
Group=fastvep
WorkingDirectory=/opt/fastvep/state

Environment=FASTVEP_PORT=3000
Environment=FASTVEP_BIND=127.0.0.1
Environment=FASTVEP_DATA_DIR=/opt/fastvep/data
# The genome that answers a request naming none; every other genome of the data directory answers when named.
Environment=FASTVEP_GFF3=/opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.115.gff3
Environment=FASTVEP_FASTA=/opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.dna.primary_assembly.fa
Environment=FASTVEP_SA_DIR=/opt/fastvep/data/hg38_ensembl_115/sa
Environment=FASTVEP_STATS_FILE=/opt/fastvep/state/stats.json
# Load the other genomes in the background at startup, so that no visitor waits for one.
Environment=FASTVEP_PRELOAD_GENOMES=true
# One annotation per core. A request is also bounded in the size of its answer (128 MiB) and its time (90 s).
Environment=FASTVEP_MAX_PARALLEL=4

ExecStart=/opt/fastvep/bin/fastvep-web
Restart=always
RestartSec=5

LimitNOFILE=65535
TasksMax=1024
# The host has 16 GB. The nine genomes take about 3.5 GiB of heap and the rest is memory-mapped files, which
# the kernel takes back under pressure. There is no MemoryHigh: the mapped pages of the databases count toward
# it, and it throttles them below what the host can hold (the rehearsal reached 11 GiB, with the old service
# beside it, and swapped 2 GiB). MemoryMax only guards the host against a runaway.
MemoryMax=13G
MemorySwapMax=2G

# The service reads its data and writes only its transcript caches and its counters.
NoNewPrivileges=true
ProtectSystem=strict
ReadWritePaths=/opt/fastvep/data /opt/fastvep/state
ProtectHome=true
PrivateTmp=true
PrivateDevices=true
ProtectKernelTunables=true
ProtectKernelModules=true
ProtectKernelLogs=true
ProtectControlGroups=true
ProtectClock=true
ProtectHostname=true
RestrictAddressFamilies=AF_INET AF_INET6 AF_UNIX
RestrictNamespaces=true
RestrictRealtime=true
RestrictSUIDSGID=true
LockPersonality=true
MemoryDenyWriteExecute=true
CapabilityBoundingSet=
AmbientCapabilities=
SystemCallArchitectures=native
SystemCallFilter=@system-service
UMask=0027

[Install]
WantedBy=multi-user.target
EOF
```

Start the service:

```bash
systemctl daemon-reload
systemctl enable fastvep
systemctl start fastvep

# Check it started
systemctl status fastvep
journalctl -u fastvep -f   # watch logs
```

What is in the unit, and what was measured about it on fastvep.org (8 November 2026):

- The sandbox (`ProtectSystem=strict` with only the data and state directories writable, no capabilities, a system call filter, no way to map memory that is writable and executable) was run before it replaced the old service: a transient unit with exactly these properties served all nine genomes on port 3001 beside the live service, with no system call refused and no restart.
- With the caches built, the default genome answers 3 s after the start and the other eight are loaded about 7 s after that; the first start on a gene model parses its GFF3 and builds the cache, which takes about 13 s for the human genome on ten cores (3.4 s for Arabidopsis release 63 on this host).
- Swapping the old service for this one took 3.6 s from the stop to the first answer.
- There is no `MemoryHigh`.
  The mapped pages of the databases count toward it, and it throttles them below what the host can hold: the rehearsal, with `MemoryHigh=11G` and the old service running beside it, reached 11 GiB within its first requests and swapped 2 GiB.
  `MemoryMax` only guards the host against a runaway.
- `FASTVEP_MAX_PARALLEL=4` is one annotation per core, and a request is bounded in the size of its answer (128 MiB) and its time (90 s): [What one request may cost](docs/API.md#what-one-request-may-cost).

For one genome only, leave out `FASTVEP_DATA_DIR` and `FASTVEP_PRELOAD_GENOMES`.

### nginx reverse proxy

nginx ends TLS, compresses the answers, and turns away a client that sends too much too fast before the server has to look at it.

The rate-limit zones and the compression types belong to the `http` context, which Ubuntu's `nginx.conf` fills from `conf.d/*.conf`.
They go in a file of their own, and not in the site file: `limit_req_zone` inside a `server` block is refused by `nginx -t` ("directive is not allowed here"; checked on nginx 1.24 of Ubuntu 24.04, which is what the earlier version of this guide had).

```bash
cat > /etc/nginx/conf.d/fastvep.conf << 'EOF'
# An annotation request takes seconds, so the rate is low and the burst is short.
limit_req_zone  $binary_remote_addr zone=fastvep_annotate:10m rate=2r/s;
limit_conn_zone $binary_remote_addr zone=fastvep_conn:10m;

# nginx compresses text/html and nothing else unless it is told: the answer of /api/annotate is JSON.
# A real answer for 200 variants was 3.25 MB, and 0.26 MB (12 times less) with gzip at nginx's default level.
gzip_types application/json text/plain text/css application/javascript image/svg+xml;
gzip_min_length 1024;
gzip_vary on;
EOF

cat > /etc/nginx/sites-available/fastvep << 'EOF'
server {
    listen 80;
    server_name your-domain.com;

    # A VCF is posted inside a JSON string, in which a tab and a newline take two bytes each, so the body is a
    # little larger than the file. The server refuses a body past 10 MiB (--max-body-size); nginx allows more, so
    # that the refusal a client sees is the server's.
    client_max_body_size 12m;

    add_header X-Content-Type-Options nosniff always;
    add_header Referrer-Policy strict-origin-when-cross-origin always;
    add_header X-Frame-Options SAMEORIGIN always;

    location = /api/annotate {
        limit_req  zone=fastvep_annotate burst=10 nodelay;
        limit_req_status 429;
        limit_conn fastvep_conn 8;
        limit_conn_status 429;

        proxy_pass http://127.0.0.1:3000;
        proxy_http_version 1.1;
        proxy_set_header Connection "";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # The server stops a request at 90 s (--max-request-seconds). nginx waits a little longer, so that the
        # client reads the server's answer and not a 504.
        proxy_read_timeout 150s;
        proxy_send_timeout 150s;
    }

    location / {
        proxy_pass http://127.0.0.1:3000;
        proxy_http_version 1.1;
        proxy_set_header Connection "";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
EOF

ln -sf /etc/nginx/sites-available/fastvep /etc/nginx/sites-enabled/
rm -f /etc/nginx/sites-enabled/default
nginx -t && systemctl reload nginx
```

This configuration was run in a container with the same nginx (1.24.0, Ubuntu 24.04): `nginx -t` passes, the 12th request of a burst of 15 to `/api/annotate` is answered `429`, a 13 MiB body is answered `413` and an 11 MiB one passes, the three headers are sent, and a JSON answer is gzipped.
Do not put an `add_header` inside a `location` of your own without repeating these three: nginx drops the inherited ones as soon as a block defines any.

### TLS with Let's Encrypt, Cloudflare, or Hostinger

**Option A: Let's Encrypt (free, auto-renewing)**

```bash
certbot --nginx -d your-domain.com
# Follow prompts; auto-renews via systemd timer
certbot renew --dry-run   # proves the renewal works before you rely on it
```

Once a dry run of the renewal has succeeded, tell browsers to use HTTPS only.
Put this line beside the other `add_header` lines, in the `server` block with `listen 443` that certbot made:

```nginx
add_header Strict-Transport-Security "max-age=31536000" always;
```

A browser that has seen it refuses plain HTTP to your domain for a year, so do not add it while the certificate can still fail.

**Option B: Hostinger (Direct DNS)**

If you bought your domain on Hostinger and want to point it directly to your Hetzner server without Cloudflare:

1. Log in to your **Hostinger Control Panel (hPanel)**.
2. Go to **Domains** and select your domain.
3. Click on **DNS / Nameservers** in the left sidebar.
4. Add or edit an **A Record**:
   - **Type**: `A`
   - **Name**: `@` (for the main domain) or `vep` (for a subdomain)
   - **Points to**: `YOUR_SERVER_IP` (your Hetzner server's public IP)
   - **TTL**: default (usually 14400)
5. Wait for propagation (usually 15-30 mins) and then run the **Let's Encrypt** setup (Option A) on your server.

**Option C: Cloudflare (recommended — free DDoS protection + CDN)**

1. **Sign up** at [cloudflare.com](https://www.cloudflare.com/) (free plan)
2. **Add your domain**: Enter your domain, Cloudflare scans existing DNS records
3. **Update nameservers**: Point your domain registrar's nameservers to Cloudflare's (shown during setup)
4. **DNS records**: Add an A record:
   - Type: `A`
   - Name: `@` (or subdomain like `vep`)
   - IPv4: your Hetzner server IP
   - Proxy: **Proxied** (orange cloud icon ON) — this enables CDN + DDoS protection
5. **SSL/TLS settings** (Cloudflare dashboard → SSL/TLS):
   - Mode: **Full (strict)**
   - Edge Certificates: enabled (automatic)
6. **Origin certificate** (SSL/TLS → Origin Server):
   - Click "Create Certificate"
   - Let Cloudflare generate key + cert (RSA 2048, 15 years)
   - Save `origin.pem` and `origin-key.pem` to server:
   ```bash
   mkdir -p /etc/ssl/cloudflare
   # paste cert into /etc/ssl/cloudflare/origin.pem
   # paste key into /etc/ssl/cloudflare/origin-key.pem
   chmod 600 /etc/ssl/cloudflare/origin-key.pem
   ```
7. **Update nginx** to use the origin certificate:
   ```nginx
   server {
       listen 443 ssl;
       server_name your-domain.com;
       ssl_certificate /etc/ssl/cloudflare/origin.pem;
       ssl_certificate_key /etc/ssl/cloudflare/origin-key.pem;

       # ... same proxy_pass config as above ...
   }
   server {
       listen 80;
       server_name your-domain.com;
       return 301 https://$host$request_uri;
   }
   ```
8. **Recommended Cloudflare settings** (free tier):
   - Caching → Cache Rules: Cache the HTML page (`/` and `/index.html`) with short TTL (1 hour)
   - Security → WAF: Managed rules enabled (automatic)
   - Speed → Auto Minify: Enable HTML/CSS/JS minification
   - Analytics: Monitor traffic and threats for free

## ACMG on the web needs the reference stack

The classifier reads nine databases ([docs/ACMG_SETUP.md](docs/ACMG_SETUP.md) lists them and what each drives).
A genome whose `sa/` directory lacks one can still be annotated, but its ACMG classes would rest on less evidence than they appear to, and nothing in a VUS says that the evidence was never looked at.
So the web server answers `"acmg": true` only for a genome that holds all nine, and refuses the request for any other with `422`, naming what is missing.
`--acmg-allow-partial` answers it anyway and says so in each answer; leave it off on a server other people use.
`/api/status` and `/api/genomes` carry the verdict (`acmg.ready`, `acmg.missing`), and the web page offers ACMG only for a genome where it will be answered.
The details are in [docs/API.md](docs/API.md#acmg-needs-the-reference-stack).

Serving ACMG for a genome is therefore a matter of building the stack into its `sa/` directory:

```bash
# Five of the nine sources: ClinVar, its protein and splice index, gnomAD gene constraints, ClinGen and
# RepeatMasker (about 1 GB of downloads; docs/SA_DATABASES.md has the manifest)
scripts/build-sa-databases.sh --assembly GRCh38 --out /opt/fastvep/data/hg38_ensembl_115/sa --sources small
# The other four: gnomAD itself, REVEL, and SpliceAI and PhyloP, which come out of the gnomAD files.
# docs/ACMG_SETUP.md builds them, and says which gnomAD release to build
python3 scripts/check_acmg_stack.py \
  --gff3 /opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.115.gff3 \
  --fasta /opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir /opt/fastvep/data/hg38_ensembl_115/sa
systemctl restart fastvep
```

`check_acmg_stack.py` annotates eight real ClinVar variants and exits non-zero when a source answers nothing, so a build that succeeded and a stack that answers are not confused.
The one cost that is large is gnomAD, which the build script keeps whole on disk: on a 160 GB disk build it one chromosome at a time and delete each VCF, as the loop of the GRCh38 section above does (the 47 GB largest genomes VCF and the 19 GB largest exomes VCF are the temporary space), or restrict it to regions with `--regions`.
GRCh37 and the other organisms have no ACMG stack here, and the server says so for each when asked.

## 5. Verify

```bash
# From your machine
curl -s https://your-domain.com/api/status | python3 -m json.tool
```

The answer is described in [docs/API.md](docs/API.md#4-verify-it-is-up).
Read four things in it before you trust an annotation:

- `transcripts` is not zero (a zero makes every variant `intergenic_variant`).
- `has_fasta` is true.
- `sa_sources` and `sa_info` list the databases you built, with the release date of each, and `sa_skipped` is empty (a file that was found and not loaded is named there).
- If you serve ACMG, `acmg.ready` is true: the server refuses `"acmg": true` for a genome that lacks a source of the reference stack, and says which.

`GET /api/genomes` says the same for each genome of the data directory.

```bash
curl -s https://your-domain.com/api/genomes | python3 -m json.tool | head -60

# Test annotation: a request names its genome ("genome": "mouse_grcm39_ensembl_115"); this one names none and gets the default
curl -s -X POST https://your-domain.com/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "17\t43124025\t.\tA\tG\t50\tPASS\t.", "pick": false}' | python3 -m json.tool
```

## 6. Maintenance

### Updates

```bash
# Build new release on dev machine
cross build --release --target x86_64-unknown-linux-gnu -p fastvep-web

# Upload beside the running binary, keep the previous one, move the new one into place, restart
scp target/x86_64-unknown-linux-gnu/release/fastvep-web root@YOUR_SERVER_IP:/opt/fastvep/bin/fastvep-web.new
ssh root@YOUR_SERVER_IP '
  set -e
  cd /opt/fastvep/bin
  chmod +x fastvep-web.new
  cp -p fastvep-web fastvep-web.previous
  mv fastvep-web.new fastvep-web
  systemctl restart fastvep
'

# Look at /api/status (version, genomes, acmg). To go back:
# ssh root@YOUR_SERVER_IP 'cd /opt/fastvep/bin && mv fastvep-web.previous fastvep-web && systemctl restart fastvep'
```

The `mv` replaces the file in one step, so a copy that is still running keeps the file it started from.
A change that you cannot take back with the previous binary, such as a new data layout, is better tried first on another port with its own `FASTVEP_STATS_FILE`, beside the running service, and swapped in once `/api/status` and a few annotations look right.

### Update SA databases

ClinVar is the one that moves.
Build the new database outside `sa/`, move it in, and restart, because a server opens its databases when it starts:

```bash
ssh root@YOUR_SERVER_IP
G=/opt/fastvep/data/hg38_ensembl_115        # and hg37_ensembl_87 with vcf_GRCh37 below
mkdir -p $G/build.tmp && cd $G/build.tmp
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh38
mv clinvar.osa2 $G/sa/clinvar.osa2          # same filesystem: the move is atomic
cd .. && rm -r build.tmp
systemctl restart fastvep
curl -s https://your-domain.com/api/status | python3 -m json.tool | grep -A3 '"clinvar"'   # the version is the release date
```

If the database you replace is in the older v1 format (`clinvar.osa` with `clinvar.osa.idx`), remove both files when the `.osa2` goes in.
A directory loads every database it holds, and when two of them know the same variant under one key the one that sorts last by path is kept ([docs/ACMG_SETUP.md](docs/ACMG_SETUP.md) says how the run reports it), so do not leave the choice to the order of two file names.

### Monitoring

```bash
# Check service health
systemctl status fastvep

# Watch logs
journalctl -u fastvep -f

# Memory usage
systemctl show fastvep -p MemoryCurrent

# Disk usage
du -sh /opt/fastvep/data/*
df -h /
```

The journal grows without anyone watching it: on fastvep.org it had reached 3.9 GB.
Give it a ceiling:

```bash
mkdir -p /etc/systemd/journald.conf.d
printf '[Journal]\nSystemMaxUse=500M\nMaxRetentionSec=1month\n' > /etc/systemd/journald.conf.d/fastvep.conf
systemctl restart systemd-journald
journalctl --disk-usage
```

### Backups

Everything that is downloaded can be downloaded again, but not as it was: ClinVar moves, and a result is only reproducible against the release it was made with.
Keep the databases of the releases you have served:

```bash
# Weekly backup of the SA databases (add to crontab)
tar czf /opt/fastvep/backup/sa_$(date +%Y%m%d).tar.gz /opt/fastvep/data/*/sa/
```

## 7. Scaling (if needed later)

| Situation | Solution | Cost |
|-----------|----------|------|
| More traffic | Upgrade to CCX33 (8 vCPU, 32GB) | ~$55/month |
| High availability | Second server + Cloudflare load balancing | ~$60/month |
| Multiple genomes | Use `--data-dir` with subdirectories per organism (nine on the 16 GB server of fastvep.org) | Same server |
| Batch WGS jobs | Run `fastvep` CLI in background, rate-limited | Same server |
| Global users | Add a second server in another region | ~$60/month |

## Architecture Diagram

```
Users (browser)
      │
      ▼
  Cloudflare (free)
  ├─ DDoS protection
  ├─ TLS termination
  └─ Caching
      │
      ▼
  Hetzner CCX23 ($30/mo)
  ├─ nginx (TLS, compression, rate limiting)
  └─ fastvep-web (port 3000)
     ├─ nine genomes, each with its own transcripts, FASTA and SA databases
     │  (GRCh38 answers a request that names none: 510k transcripts, Ensembl 115)
     ├─ FASTA references (memory-mapped)
     └─ SA databases of GRCh38, as many as you built
        ├─ ClinVar
        ├─ gnomAD v4
        ├─ dbSNP
        ├─ 1000 Genomes
        ├─ REVEL
        └─ SpliceAI
```

#!/usr/bin/env bash

# ==============================================================================
# SCRIPT: Preparation of Brassica napus Darmor-bzh v10 reference proteome
#         and bidirectional BLASTp analysis
#
# Description:
# This script prepares the Brassica napus cv. Darmor-bzh v10 reference proteome
# for identification of homoeologous gene pairs between the BnA and BnC
# subgenomes.
#
# The workflow:
#   1. Decompresses the reference protein FASTA file.
#   2. Separates protein sequences belonging to the BnA and BnC subgenomes.
#   3. Retains only the first protein isoform (.1).
#   4. Creates BLAST protein databases for each subgenome.
#   5. Performs bidirectional BLASTp searches (BnA vs BnC and BnC vs BnA).
#
# Reference genome:
# Brassica napus cv. Darmor-bzh v10
#
# Required reference file:
# BnapusDarmor-bzh_proteins.fasta.gz
#
# Reference data source:
# BnaOmics
# https://bnaomics.ocri-genomics.net/download/public/
# reference-genome-sequences-and-gene-annotation/
# B_napus_cv_Darmor_bzh/v10/
# ==============================================================================


# ==============================================================================
# 0. Requirements
# ==============================================================================

# Required command-line tools:
#   - BLAST+ (blastp and makeblastdb)
#   - pv (used to monitor progress during BLAST searches)
#
# These tools must be installed before running the script.
# Installation instructions are provided in the repository README.

command -v blastp >/dev/null 2>&1 || {
    echo "Error: blastp was not found in PATH." >&2
    exit 1
}

command -v makeblastdb >/dev/null 2>&1 || {
    echo "Error: makeblastdb was not found in PATH." >&2
    exit 1
}

command -v pv >/dev/null 2>&1 || {
    echo "Error: pv was not found in PATH." >&2
    exit 1
}

# Display installed BLASTp version.
blastp -version


# ==============================================================================
# 1. Define input and output directories
# ==============================================================================

# REFERENCE_DIR:
# Directory containing the compressed Darmor-bzh v10 protein FASTA file.
#
# WORK_DIR:
# Directory where intermediate files, BLAST databases and BLAST results
# will be generated.
#
# NUM_THREADS:
# Number of processor threads used by BLASTp.
#
# These values can be modified here or supplied as environment variables
# when running the script.

REFERENCE_DIR="${REFERENCE_DIR:-./reference}"
WORK_DIR="${WORK_DIR:-./Homologous_Analysis_DarmorV10}"
NUM_THREADS="${NUM_THREADS:-8}"


# Check that the reference directory exists.
if [ ! -d "$REFERENCE_DIR" ]; then
    echo "Error: reference directory not found: $REFERENCE_DIR" >&2
    exit 1
fi


# Convert the reference directory to an absolute path.
REFERENCE_DIR="$(cd "$REFERENCE_DIR" && pwd)"


# Define the reference protein FASTA file.
PROTEIN_FASTA="${REFERENCE_DIR}/BnapusDarmor-bzh_proteins.fasta.gz"


# Check that the reference protein FASTA file exists.
if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "Error: reference protein FASTA file not found:" >&2
    echo "$PROTEIN_FASTA" >&2
    exit 1
fi


# Create the directory structure used throughout the analysis.
mkdir -p "${WORK_DIR}"/{proteins,db,results,tmp}


# Convert the working directory to an absolute path.
WORK_DIR="$(cd "$WORK_DIR" && pwd)"


# ==============================================================================
# 2. Decompress the reference protein FASTA file
# ==============================================================================

cd "$WORK_DIR" || exit 1

gunzip -c "$PROTEIN_FASTA" > proteins/all_proteins.fa


# ==============================================================================
# 3. Separate BnA and BnC subgenomes
# ==============================================================================

# Protein sequences whose FASTA headers begin with "A" are assigned to
# the BnA subgenome.

awk '/^>/ {keep = ($0 ~ /^>A/)} keep {print}' \
    proteins/all_proteins.fa > proteins/A.fa


# Protein sequences whose FASTA headers begin with "C" are assigned to
# the BnC subgenome.

awk '/^>/ {keep = ($0 ~ /^>C/)} keep {print}' \
    proteins/all_proteins.fa > proteins/C.fa


# ==============================================================================
# 3.1. Sanity check after subgenome separation
# ==============================================================================

# Count the number of protein sequences assigned to each subgenome.

echo "Number of BnA protein sequences:"
grep -c "^>" proteins/A.fa

echo "Number of BnC protein sequences:"
grep -c "^>" proteins/C.fa


# ==============================================================================
# 4. Retain only the first protein isoform (.1)
# ==============================================================================

# To avoid redundancy caused by alternative protein isoforms, only sequences
# corresponding to the first isoform (.1) are retained.

for genome in A C; do
    awk '/^>/ {keep = ($0 ~ /\.1/)} keep {print}' \
        "proteins/${genome}.fa" > "proteins/${genome}_iso1.fa"
done


# ==============================================================================
# 4.1. Sanity check after isoform filtering
# ==============================================================================

# Count the number of protein sequences retained after filtering.

echo "Number of BnA protein sequences after isoform filtering:"
grep -c "^>" proteins/A_iso1.fa

echo "Number of BnC protein sequences after isoform filtering:"
grep -c "^>" proteins/C_iso1.fa


# ==============================================================================
# 5. Create BLAST protein databases
# ==============================================================================

# Create a BLAST protein database for the BnA subgenome.

makeblastdb \
    -in proteins/A_iso1.fa \
    -dbtype prot \
    -out db/A


# Create a BLAST protein database for the BnC subgenome.

makeblastdb \
    -in proteins/C_iso1.fa \
    -dbtype prot \
    -out db/C


# ==============================================================================
# 6. Run bidirectional BLASTp searches
# ==============================================================================

# BnA proteins queried against the BnC protein database.

pv proteins/A_iso1.fa | \
blastp \
    -query - \
    -db db/C \
    -out results/A_vs_C.tsv \
    -outfmt 6 \
    -num_threads "$NUM_THREADS"


# BnC proteins queried against the BnA protein database.

pv proteins/C_iso1.fa | \
blastp \
    -query - \
    -db db/A \
    -out results/C_vs_A.tsv \
    -outfmt 6 \
    -num_threads "$NUM_THREADS"
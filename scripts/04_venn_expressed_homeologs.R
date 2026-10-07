# ==============================================================================
# SCRIPT: Venn diagrams of expressed BnA and BnC genes based on TPM
#
# Description:
# This script identifies expressed genes from the BnA and BnC subgenomes of
# Brassica napus and generates separate Venn diagrams showing their overlap
# across the defined lineages.
#
# A gene is considered expressed within a lineage when TPM >= 1 in at least one
# sample belonging to that lineage.
#
# IMPORTANT:
# TPM sample columns must begin with the corresponding lineage prefix because
# the script selects lineage-specific samples using starts_with(lineage).
#
# For example, if a lineage is defined as "Line1", valid sample names include:
#   Line1_rep1
#   Line1_rep2
#   Line1_UP_rep1
#   Line1_P_rep1
#
# Names such as "rep1_Line1" will not be selected by the current workflow.
#
# The first column of the TPM matrix must contain gene identifiers. The script
# renames this column to "gene.id" internally.
# ==============================================================================


# ==============================================================================
# 1. Load required packages
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
})


# ==============================================================================
# 2. Define input and output paths
# ==============================================================================

# TPM_FILE:
# TPM-normalized expression matrix in Excel format.
#
# OUTPUT_DIR:
# Directory where the Venn diagrams will be saved.
#
# These paths can be modified here or supplied as environment variables when
# running the script.

TPM_FILE <- Sys.getenv(
  "TPM_FILE",
  unset = "./data/tpm_counts.xlsx"
)

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  unset = "./results/Venn_expressed_genes"
)


# Check that the TPM file exists.
if (!file.exists(TPM_FILE)) {
  stop(
    paste0(
      "TPM file not found: ",
      TPM_FILE,
      "\nSet TPM_FILE to the correct path before running the script."
    )
  )
}


# Create output directory.
dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================================
# 3. Load TPM data
# ==============================================================================

message(">>> Reading TPM matrix...")

tpm <- read_xlsx(TPM_FILE)

# The first column is assumed to contain gene identifiers.
colnames(tpm)[1] <- "gene.id"


# ==============================================================================
# 4. Define lineages
# ==============================================================================

# Replace these example lineage names with the prefixes used in the TPM matrix.
#
# Each lineage name must match the beginning of the corresponding TPM sample
# columns because lineage-specific samples are selected using starts_with().

lineages <- c("Line1", "Line2", "Line3")


# ==============================================================================
# 5. Separate BnA and BnC subgenomes
# ==============================================================================

# Gene identifiers beginning with "A" are assigned to the BnA subgenome.
tpm_A <- tpm %>%
  filter(grepl("^A", gene.id))

# Gene identifiers beginning with "C" are assigned to the BnC subgenome.
tpm_C <- tpm %>%
  filter(grepl("^C", gene.id))


# ==============================================================================
# 6. Identify expressed genes within each lineage
# ==============================================================================

# A gene is retained when TPM >= 1 in at least one sample belonging to the
# corresponding lineage.

get_expressed_genes <- function(df, lineage) {
  
  lineage_data <- df %>%
    select(gene.id, starts_with(lineage))
  
  # Stop with an informative message if no TPM columns match the lineage prefix.
  if (ncol(lineage_data) == 1) {
    stop(
      paste0(
        "No TPM sample columns were found for lineage '",
        lineage,
        "'. Sample column names must begin with the lineage prefix."
      )
    )
  }
  
  lineage_data %>%
    filter(rowSums(select(., -gene.id) >= 1) > 0) %>%
    pull(gene.id)
}


# Build one list of expressed genes for BnA and one for BnC.
A_list <- setNames(
  lapply(
    lineages,
    function(lineage) get_expressed_genes(tpm_A, lineage)
  ),
  lineages
)

C_list <- setNames(
  lapply(
    lineages,
    function(lineage) get_expressed_genes(tpm_C, lineage)
  ),
  lineages
)


# ==============================================================================
# 7. Generate Venn diagram for the BnA subgenome
# ==============================================================================

p_A <- ggVennDiagram(
  A_list,
  label_alpha = 0,
  set_color = c("#104E8B", "#fee08b", "#CD0000"),
  edge_size = 0.5
) +
  scale_fill_gradient(
    low = "#E3F2FD",
    high = "#00BFFF"
  ) +
  ggtitle("Genes Subgenome A - TPM >= 1") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Venn_Subgenome_A_TPM.png"
  ),
  plot = p_A,
  width = 6,
  height = 6
)


# ==============================================================================
# 8. Generate Venn diagram for the BnC subgenome
# ==============================================================================

p_C <- ggVennDiagram(
  C_list,
  label_alpha = 0,
  set_color = c("#104E8B", "#fee08b", "#CD0000"),
  edge_size = 0.5
) +
  scale_fill_gradient(
    low = "#FFF8DC",
    high = "#9BCD9B"
  ) +
  ggtitle("Genes Subgenome C - TPM >= 1") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Venn_Subgenome_C_TPM.png"
  ),
  plot = p_C,
  width = 6,
  height = 6
)

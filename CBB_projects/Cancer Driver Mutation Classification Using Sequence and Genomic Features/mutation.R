library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query, files.per.chunk = 10)

library(maftools)

library(dplyr)

# Function to read and coerce types to character
safe_read_maf <- function(f) {
  tryCatch({
    df <- read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
    df[] <- lapply(df, as.character)  # convert all columns to character
    return(df)
  }, error = function(e) NULL)
}

# Apply across all MAF files
all_mafs <- bind_rows(lapply(maf_paths, safe_read_maf))


selected <- all_mafs[, c("Hugo_Symbol", "Chromosome", "Start_Position",
                         "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification")]

colnames(selected) <- c("Gene", "Chromosome", "Position", "Ref", "Alt", "Classification")

write.csv(selected, "tcga_brca_mutations.csv", row.names = FALSE)

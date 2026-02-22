########################################################
#SIEVE: One-stop differential expression, variability,
#       and skewness analyses using RNA-Seq data
#Authors: Hongxiang Li and Dr. Tsung Fei Khang
#Email: hxli@ynnu.edu.cn (H. Li)
#Date: 22 Feb. 2026
#Part 9: R Codes for cross-methodologcial and cross-data
#validation: comparison of ORA and GSEA
########################################################

## R packages downloaded
#devtools::install_github("alserglab/fgsea")
#BiocManager::install("org.Hs.eg.db")
library(fgsea)
library(org.Hs.eg.db)

###########
# Read data
# GSE249477_raw_count_normalize_04-10-2025.csv.gz
# Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249477
# Nakayama data
Data_AD = read.csv(gzfile('GSE249477_raw_count_normalize_04-10-2025.csv.gz'),
                   header = T, check.names = F, row.names = NULL)

# Extract Total counts columns
count_cols <- grep("Total counts", colnames(Data_AD), value = TRUE)

# Specify AD and Control samples
# AD: DK22119_1 TO DK22119_21
# Control: DK22119_43 TO DK22119_63
ad_samples <- paste0("DK22119_", sprintf("%02d", 1:21))
ctrl_samples <- paste0("DK22119_", sprintf("%02d", 43:63))
selected_counts <- count_cols[
  grepl(paste(c(ctrl_samples, ad_samples), collapse="|"), count_cols)
]

# Extract required columns
Data_subset <- Data_AD[, c("Identifier", "Feature ID", selected_counts)]

# Keep valid Ensembl IDs only
Data_subset <- Data_subset[
  grepl("^ENSG", Data_subset$Identifier),
]

# Remove duplicated Identifiers
Data_subset <- Data_subset[
  !duplicated(Data_subset$Identifier),
]

# Clean column names
colnames(Data_subset)[-(1:2)] <- gsub(
  " \\(GE\\) - Total counts",
  "",
  colnames(Data_subset)[-(1:2)]
)

# Set row names
rownames(Data_subset) <- Data_subset$Identifier
Data_subset <- Data_subset[, -(1:2)]
dim(Data_subset) # 21479    42
Data_subset <- na.omit(Data_subset)
dim(Data_subset) # 21466    42
Data_subset <- cbind(Data_subset[, 22:42], Data_subset[, 1:21])


###########
# Filtering
###########
CPM <- cpm(Data_subset)
keep <- rowMeans(CPM[, 1:21]) > 0.5 &
  rowMeans(CPM[,22:42]) > 0.5 &
  apply(Data_subset[, 1:21], 1, function(k) length(k[k == 0])/length(k)) < 0.85 &
  apply(Data_subset[, 22:42], 1, function(k) length(k[k == 0])/length(k)) < 0.85
AD_control_filter <- Data_subset[keep, ]
# dim(AD_control_filter) # 13078  42


#############################
### DE/DV/DS Test using SIEVE
#############################
# clr transformation
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}

## DE/DV/DS test
t1 <- proc.time()
clr_counts <- clr.transform(data = AD_control_filter)
group2 = c(rep(0, 21), rep(1, 21)) # 21 controls vs 21 AD
clrSeq_result <- clrSeq(clr_counts, group = group2)
clrSIEVE_result <- clrSIEVE(clrSeq_result = clrSeq_result,
                            alpha_level = 0.05,
                            order_DE = F,
                            order_LFC = F,
                            order_DS = F,
                            order_sieve = F)
as.numeric(proc.time() - t1)[3]
# around 25 minutes

########
### GSEA
########

# Extract DE/DV/DS test results
DE_AD_table <- clrSIEVE_result$clrDE_test
DV_AD_table <- clrSIEVE_result$clrDV_test
DS_AD_table <- clrSIEVE_result$clrDS_test

# Check NA values
sum(is.na(DE_AD_table))
# DE_AD_table[apply(DE_AD_table, 1, function(x) any(is.na(x))), ]
sum(is.na(DV_AD_table))
# DV_AD_table[apply(DV_AD_table, 1, function(x) any(is.na(x))), ]
sum(is.na(DS_AD_table))
# DS_AD_table[apply(DS_AD_table, 1, function(x) any(is.na(x))), ]

DE_AD_table <- na.omit(DE_AD_table)
DV_AD_table <- na.omit(DV_AD_table)
DS_AD_table <- na.omit(DS_AD_table)

# Calculate composite ranking score
DE_DV_DS_Wald_Stat_Mat <- cbind(DE_AD_table$DE,
                                DV_AD_table$LFC,  # for DV
                                DS_AD_table$DS)

# range(DE_AD_table$DE)
# range(DV_AD_table$LFC)
# range(DS_AD_table$DS)

row.names(DE_DV_DS_Wald_Stat_Mat) <- rownames(DE_AD_table)

# Identify index of maximum absolute statistic per gene
idx <- apply(abs(DE_DV_DS_Wald_Stat_Mat), 1, which.max)

# Extract corresponding statistic
max_vec <- DE_DV_DS_Wald_Stat_Mat[
  cbind(seq_len(nrow(DE_DV_DS_Wald_Stat_Mat)), idx)
]
names(max_vec) <- rownames(DE_DV_DS_Wald_Stat_Mat)

max_vec_ranked <- sort(max_vec, decreasing = T)

max_vec_ranked <- max_vec_ranked[
  names(max_vec_ranked) %in%
    names(which(table(names(max_vec_ranked)) == 1))
]

head(max_vec_ranked)
tail(max_vec_ranked)

# Reference list of 18 GO IDs
go_ids_18 <- c(
  "GO:0034613", "GO:0072599", "GO:0090150", "GO:0045047",
  "GO:0006614", "GO:0022610", "GO:0031589", "GO:0010647",
  "GO:0070848", "GO:0030198", "GO:0001568", "GO:0042060",
  "GO:0006412", "GO:0006401", "GO:0019884", "GO:0045619",
  "GO:0019363", "GO:0019752"
)

# GO descriptions for Supplementary Table
go_descriptions <- c(
  "GO:0034613" = "cellular protein localization",
  "GO:0072599" = "establishment of protein localization to endoplasmic reticulum",
  "GO:0090150" = "establishment of protein localization to membrane",
  "GO:0045047" = "protein targeting to ER",
  "GO:0006614" = "SRP-dependent cotranslational protein targeting to membrane",
  "GO:0022610" = "biological adhesion",
  "GO:0031589" = "cell-substrate adhesion",
  "GO:0010647" = "positive regulation of cell communication",
  "GO:0070848" = "response to growth factor",
  "GO:0030198" = "extracellular matrix organization",
  "GO:0001568" = "blood vessel development",
  "GO:0042060" = "wound healing",
  "GO:0006412" = "translation",
  "GO:0006401" = "RNA catabolic process",
  "GO:0019884" = "antigen processing and presentation of exogenous antigen",
  "GO:0045619" = "regulation of lymphocyte differentiation",
  "GO:0019363" = "pyridine nucleotide biosynthetic process",
  "GO:0019752" = "carboxylic acid metabolic process"
)

###########################################
### Build Gene Sets (Mapping GO to ENSEMBL)
###########################################

# Primary mapping: GO to ENSEMBL
go_map_main <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_ids_18,
  keytype = "GOALL",
  columns = c("GOALL", "ENSEMBL")
) %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::distinct()

# Split into named list for fgsea
gene_sets_table <- split(go_map_main$ENSEMBL, go_map_main$GOALL)

# Remedial mapping for missing GO IDs
missing_ids <- setdiff(go_ids_18, names(gene_sets_table))
if (length(missing_ids) > 0) {
  for (go_id in missing_ids) {
    entrez_vec <- tryCatch({
      obj <- mget(go_id, envir = org.Hs.egGO2ALLEGS, ifnotfound = NA)[[1]]
      if (all(is.na(obj))) character(0) else as.character(obj)
    }, error = function(e) character(0))

    if (length(entrez_vec) == 0) next

    ens_map <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = entrez_vec,
      column = "ENSEMBL", keytype = "ENTREZID", multiVals = "list"
    )

    ens_vec <- unique(unlist(ens_map))
    ens_vec <- ens_vec[!is.na(ens_vec)]

    if (length(ens_vec) > 0) {
      gene_sets_table[[go_id]] <- ens_vec
    }
  }
}

###########################################
### GSEA calculation and results formatting
###########################################

set.seed(100)
GSEA_DE_DV_DS <- fgseaSimple(
  pathways = gene_sets_table,
  stats = max_vec_ranked,
  nperm = 5*10^5,
  minSize = 5,
  maxSize = 2000,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL
)

# Format final data frame for Table 5
final_output <- as.data.frame(GSEA_DE_DV_DS) %>%
  dplyr::filter(pathway %in% names(go_descriptions)) %>%
  dplyr::mutate(
    pathway = as.character(pathway),
    description = go_descriptions[pathway]
  ) %>%
  dplyr::select(
    GO_ID = pathway,
    description,
    NES,
    padj
  ) %>%
  dplyr::arrange(padj)

# Display final table
print(final_output)

### END ###

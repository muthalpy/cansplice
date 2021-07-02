library(tidyverse)
library(SummarizedExperiment)
library(edgeR)

dir_gdc <- file.path(getwd(), "Data","Queried")
dir_mda <- file.path(getwd(), "Data", "TCGAspliceseqData", "spliceseqdata")
dir_app <- file.path(getwd(), "App",  "datasets")

sites <- c("BRCA", "CESC", "COAD", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "UCEC")

for (site in sites){

load(file.path(dir_gdc, paste0("TCGA-", site, "-queried-expression.RData")))
rse_expr <- eval(parse(text = paste0("rse_", site)))

pdat_gdc <- as.data.frame(colData(rse_expr)) %>%
  rownames_to_column(var = "barcode2") %>%
  transmute(barcode2,
            patient= ifelse(definition == "Solid Tissue Normal", paste0(patient, "-", "Norm"), patient ), 
            barcode, 
            definition,
            depth = colSums(assay(rse_expr)[, barcode2]),
            race,
            status = factor(vital_status,
                            levels = c("Alive", "Dead"),
                            labels = c("0", "1")),
            month = ifelse(status == "0", 
                           days_to_last_follow_up/365.25*12, 
                           days_to_death/365.25*12),
            subtype_BRCA_Subtype_PAM50=  ifelse(site == "BRCA",factor(subtype_BRCA_Subtype_PAM50), ""),
            subtype_PRAD_Subtype_Gleason=  ifelse(site == "PRAD",factor(subtype_BRCA_Subtype_PAM50), "")
            ) %>%
  group_by(patient) %>%
  filter(depth == max(depth)) %>%
  ungroup()


###################
rse_expr <- rse_expr[, !is.na(rse_expr$race)]
adat_count <- assay(rse_expr)
design <- model.matrix(~ definition + race, 
                       data = as.data.frame(colData(rse_expr)))
norm_factor <- calcNormFactors(adat_count)
lib_size <- colSums(adat_count) * norm_factor
adat_expr <- voom(adat_count, 
                  lib.size = lib_size, 
                  design = design)$E

rse_expr <- SummarizedExperiment(assays = adat_expr,
                                 colData = colData(rse_expr),
                                 rowData = rowData(rse_expr))
#########################

 pdat <- pdat_gdc

file_psi <- file.path(dir_mda, paste0("PSI_download_", site, ".zip"))
psi <- read.table(unzip(file_psi), 
                  header = T, 
                  na.strings = "null", 
                  stringsAsFactors = F) %>%
  mutate(dse_id = paste0(gsub("-", "", symbol), ".", as_id))

cdat_psi <- pdat %>%
  column_to_rownames(var = "patient")

adat_psi <- psi %>%
  column_to_rownames(var = "dse_id") %>%
  purrr::set_names(gsub("_", "-", colnames(.))) %>%
  dplyr::select(one_of(rownames(cdat_psi))) %>%
  as.matrix() 
rdat_psi <- psi %>%
  mutate(pct_with_values_subcols = rowMeans(!is.na(adat_psi))) %>%
  dplyr::select(-starts_with("TCGA")) %>%
  column_to_rownames(var = "dse_id") 

cids <- intersect(colnames(adat_psi), rownames(cdat_psi))
rids <- intersect(rownames(adat_psi), rownames(rdat_psi))
rse_psi <- SummarizedExperiment(assays = adat_psi[rids, cids],
                                colData = cdat_psi[cids, ],
                                rowData = rdat_psi[rids, ])

cids_cm <- intersect(colData(rse_expr)$barcode, colData(rse_psi)$barcode)
colnames(rse_expr) <- colData(rse_expr)$barcode
colnames(rse_psi) <- colData(rse_psi)$barcode
rse_expr <- rse_expr[, cids_cm]
rse_psi <- rse_psi[, cids_cm]
ncol(rse_expr)
ncol(rse_psi)

assign(paste0("rse_psi_", site), rse_psi)
save(list = paste0("rse_psi_", site),
     file = file.path(dir_app, paste0("TCGA-", site, "-psi.RData")))
rm(rse_psi)
rm(list = paste0("rse_psi_", site))
rm(adat_psi)
assign(paste0("rse_expr_", site), rse_expr)
save(list = paste0("rse_expr_", site), 
     file = file.path(dir_app, paste0("TCGA-", site, "-expression.RData")))
rm(rse_expr)
rm(list = paste0("rse_expr_", site))
rm(adat_expr)
rm(list = paste0("rse_", site))
}

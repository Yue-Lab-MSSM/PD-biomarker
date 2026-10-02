# =========================================================
# PPMI -> LCC external validation (URINE)
# Main task: HC vs LRRK2 PD
# Both files are CSV
# Both datasets already log2 transformed
# Models:
#   1) Age + sex only
#   2) Elastic net proteins only
#   3) Elastic net proteins + age + sex
# =========================================================

library(glmnet)
library(pROC)
library(ggplot2)
library(yardstick)

set.seed(23)

# =========================
# 1. FILE PATHS
# =========================
ppmi_file <- "../All_modified_log_imputation.csv"
lcc_file  <- "../The urine proteomic data of LCC.csv"

outdir <- "./"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2. HELPER FUNCTIONS
# =========================

normalize_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("m", "male", "1")] <- "male"
  x[x %in% c("f", "female", "0")] <- "female"
  factor(x, levels = c("female", "male"))
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- trimws(x)
  x
}

clean_upper <- function(x) {
  toupper(clean_text(x))
}

# treat as numeric
sanitize_numeric_cols <- function(df, cols) {
  df2 <- df
  cols <- cols[!is.na(cols) & cols != ""]
  cols <- intersect(cols, colnames(df2))
  
  for (cc in cols) {
    x <- trimws(as.character(df2[[cc]]))
    x[x %in% c("", "NA", "N/A", "NaN", "NULL")] <- NA
    x <- gsub(",", "", x)
    df2[[cc]] <- suppressWarnings(as.numeric(x))
  }
  df2
}

# row-wise median normalization
median_normalize_rows <- function(df, prot_cols) {
  df2 <- df
  mat <- as.matrix(df2[, prot_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  row_meds <- apply(mat, 1, median, na.rm = TRUE)
  mat_norm <- sweep(mat, 1, row_meds, FUN = "-")
  df2[, prot_cols] <- mat_norm
  df2
}

# impute numeric columns using TRAINING medians
median_impute_by_train <- function(train_df, test_df, cols) {
  train2 <- train_df
  test2  <- test_df
  
  for (cc in cols) {
    med <- median(train2[[cc]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    train2[[cc]][is.na(train2[[cc]])] <- med
    test2[[cc]][is.na(test2[[cc]])]   <- med
  }
  list(train = train2, test = test2)
}

# z-score using TRAINING mean/sd
scale_by_train <- function(train_df, test_df, cols) {
  train2 <- train_df
  test2  <- test_df
  
  for (cc in cols) {
    mu <- mean(train2[[cc]], na.rm = TRUE)
    sdv <- sd(train2[[cc]], na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) sdv <- 1
    train2[[cc]] <- (train2[[cc]] - mu) / sdv
    test2[[cc]]  <- (test2[[cc]] - mu) / sdv
  }
  list(train = train2, test = test2)
}

calc_auc_ci <- function(truth, prob, negative_class = "HC", positive_class = "LRRK2_PD") {
  truth <- factor(as.character(truth), levels = c(negative_class, positive_class))
  prob  <- as.numeric(prob)
  
  if (sum(truth == negative_class, na.rm = TRUE) == 0) {
    stop(paste0("No control samples found: ", negative_class))
  }
  if (sum(truth == positive_class, na.rm = TRUE) == 0) {
    stop(paste0("No case samples found: ", positive_class))
  }
  
  roc_obj <- pROC::roc(
    response = truth,
    predictor = prob,
    levels = c(negative_class, positive_class),
    direction = "<"
  )
  ci_obj <- pROC::ci.auc(roc_obj)
  
  data.frame(
    AUC = as.numeric(pROC::auc(roc_obj)),
    CI_lower = as.numeric(ci_obj[1]),
    CI_mid   = as.numeric(ci_obj[2]),
    CI_upper = as.numeric(ci_obj[3])
  )
}

plot_roc_save <- function(truth, prob, file, title_text,
                          negative_class = "HC", positive_class = "LRRK2_PD") {
  truth <- factor(as.character(truth), levels = c(negative_class, positive_class))
  prob  <- as.numeric(prob)
  
  roc_obj <- pROC::roc(
    response = truth,
    predictor = prob,
    levels = c(negative_class, positive_class),
    direction = "<"
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  x <- 1 - roc_obj$specificities
  y <- roc_obj$sensitivities
  
  pdf(file, width = 5, height = 5)
  plot(
    x, y,
    type = "s",
    xlim = c(0, 1),
    ylim = c(0, 1),
    col = "#d33524",
    lwd = 2,
    xlab = "1 - Specificity",
    ylab = "Sensitivity",
    main = title_text
  )
  abline(a = 0, b = 1, lty = 2, col = "grey")
  text(
    x = 0.72, y = 0.06,
    labels = paste0("AUC: ", sprintf("%.3f", auc_val)),
    col = "#d33524",
    cex = 1.2,
    adj = c(0, 0)
  )
  dev.off()
}

# urine protein-name cleaning
strip_protein_name <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^P151_CSF_", "", x)
  x <- sub("^P151_URINE_", "", x)
  
  x <- strsplit(x, "[;|/,]")[[1]][1]
  x <- trimws(x)
  
  x <- strsplit(x, "_")[[1]][1]
  x <- trimws(x)
  
  if (is.na(x) || x == "") return(NA_character_)
  x
}

# =========================
# 3. LOAD DATA
# =========================
ppmi <- read.csv(ppmi_file, check.names = FALSE, stringsAsFactors = FALSE)
lcc  <- read.csv(lcc_file,  check.names = FALSE, stringsAsFactors = FALSE)

# =========================
# 4. DEFINE GROUPS / COVARIATES
# =========================

# ---- PPMI ----
ppmi$group2 <- NA_character_
ppmi$group2[trimws(as.character(ppmi$GROUP)) == "HC"] <- "HC"
ppmi$group2[trimws(as.character(ppmi$GROUP)) == "PD_LRRK2"] <- "LRRK2_PD"

ppmi <- ppmi[!is.na(ppmi$group2), , drop = FALSE]
ppmi$group2 <- factor(ppmi$group2, levels = c("HC", "LRRK2_PD"))

ppmi$sex2 <- factor(
  ifelse(as.numeric(ppmi$SEX) == 1, "male", "female"),
  levels = c("female", "male")
)
ppmi$age2 <- as.numeric(ppmi$AGE)

# ---- LCC ----
cat("Raw LCC Group counts:\n")
print(table(lcc$Group, useNA = "ifany"))

lcc$group2 <- NA_character_
lcc$group2[trimws(as.character(lcc$Group)) == "HC"] <- "HC"
lcc$group2[trimws(as.character(lcc$Group)) == "LRRK2_PD"] <- "LRRK2_PD"

lcc <- lcc[!is.na(lcc$group2), , drop = FALSE]
lcc$group2 <- factor(lcc$group2, levels = c("HC", "LRRK2_PD"))

lcc$sex2 <- normalize_sex(lcc$Gender)
lcc$age2 <- as.numeric(lcc$`Age at sampe collection`)

cat("PPMI counts:\n")
print(table(ppmi$group2, useNA = "ifany"))
cat("LCC counts after filtering:\n")
print(table(lcc$group2, useNA = "ifany"))
print(table(lcc$sex2, useNA = "ifany"))

# =========================
# 5. DEFINE PROTEIN COLUMNS
# =========================

# PPMI urine: first protein column is A0A0G2JRQ6
ppmi_start <- match("A0A0G2JRQ6", colnames(ppmi))
if (is.na(ppmi_start)) stop("Could not find 'A0A0G2JRQ6' in PPMI columns.")

ppmi_meta_added_cols <- c("group2", "sex2", "age2")
ppmi_prot_cols_raw <- setdiff(colnames(ppmi)[ppmi_start:ncol(ppmi)], ppmi_meta_added_cols)

# LCC urine: first protein column is IGLV4-69
lcc_start <- match("IGLV4-69", colnames(lcc))
if (is.na(lcc_start)) stop("Could not find 'IGLV4-69' in LCC columns.")

meta_added_cols <- c("group2", "sex2", "age2")
lcc_prot_cols_raw <- setdiff(colnames(lcc)[lcc_start:ncol(lcc)], meta_added_cols)

# clean names
ppmi_gene_names <- sapply(ppmi_prot_cols_raw, strip_protein_name)
lcc_gene_names  <- sapply(lcc_prot_cols_raw,  strip_protein_name)

# make cleaned duplicate names unique
ppmi_gene_names <- make.unique(ppmi_gene_names, sep = "_dup")
lcc_gene_names  <- make.unique(lcc_gene_names,  sep = "_dup")

# write cleaned names back
colnames(ppmi)[match(ppmi_prot_cols_raw, colnames(ppmi))] <- ppmi_gene_names
colnames(lcc)[match(lcc_prot_cols_raw, colnames(lcc))]   <- lcc_gene_names

# reset protein column vectors after renaming
ppmi_prot_cols <- ppmi_gene_names
lcc_prot_cols  <- lcc_gene_names

cat("Does PPMI protein list still contain helper cols?\n")
print(intersect(ppmi_prot_cols, c("group2", "sex2", "age2")))

cat("Does LCC protein list still contain helper cols?\n")
print(intersect(lcc_prot_cols, c("group2", "sex2", "age2")))

cat("PPMI protein names not found in df after renaming:\n")
print(setdiff(ppmi_prot_cols, colnames(ppmi)))

cat("LCC protein names not found in df after renaming:\n")
print(setdiff(lcc_prot_cols, colnames(lcc)))

cat("Blank/NA PPMI protein names:\n")
print(ppmi_prot_cols[is.na(ppmi_prot_cols) | ppmi_prot_cols == ""])

cat("Blank/NA LCC protein names:\n")
print(lcc_prot_cols[is.na(lcc_prot_cols) | lcc_prot_cols == ""])

cat("Duplicated cleaned LCC protein names renamed as:\n")
print(lcc_prot_cols[grepl("_dup", lcc_prot_cols)])

cat("First 20 cleaned PPMI protein names:\n")
print(head(ppmi_prot_cols, 20))
cat("First 20 cleaned LCC protein names:\n")
print(head(lcc_prot_cols, 20))

ppmi_prot_cols <- ppmi_prot_cols[!is.na(ppmi_prot_cols) & ppmi_prot_cols != ""]
lcc_prot_cols  <- lcc_prot_cols[!is.na(lcc_prot_cols) & lcc_prot_cols != ""]

ppmi_prot_cols <- intersect(ppmi_prot_cols, colnames(ppmi))
lcc_prot_cols  <- intersect(lcc_prot_cols, colnames(lcc))

# =========================
# 6. PREPROCESS
# =========================

# both datasets already log2-transformed
ppmi <- sanitize_numeric_cols(ppmi, ppmi_prot_cols)
lcc  <- sanitize_numeric_cols(lcc,  lcc_prot_cols)

# row-wise median normalization on BOTH cohorts
ppmi <- median_normalize_rows(ppmi, ppmi_prot_cols)
lcc  <- median_normalize_rows(lcc,  lcc_prot_cols)

cat("Any all-NA PPMI protein columns after preprocessing?\n")
all_na_cols_ppmi <- names(which(colSums(!is.na(ppmi[, ppmi_prot_cols, drop = FALSE])) == 0))
print(length(all_na_cols_ppmi))
print(all_na_cols_ppmi)

cat("Any all-NA LCC protein columns after preprocessing?\n")
all_na_cols <- names(which(colSums(!is.na(lcc[, lcc_prot_cols, drop = FALSE])) == 0))
print(length(all_na_cols))
print(all_na_cols)

# =========================
# 7. HARMONIZE COMMON PROTEINS
# =========================
common_genes <- intersect(ppmi_prot_cols, lcc_prot_cols)
common_genes <- sort(common_genes)

cat("Overlapping proteins between PPMI and LCC:", length(common_genes), "\n")

ppmi_h <- data.frame(
  group2 = ppmi$group2,
  age2   = ppmi$age2,
  sex2   = ppmi$sex2,
  ppmi[, common_genes, drop = FALSE],
  check.names = FALSE
)

lcc_h <- data.frame(
  group2 = lcc$group2,
  age2   = lcc$age2,
  sex2   = lcc$sex2,
  lcc[, common_genes, drop = FALSE],
  check.names = FALSE
)

ppmi_h$group2 <- factor(as.character(ppmi_h$group2), levels = c("HC", "LRRK2_PD"))
lcc_h$group2  <- factor(as.character(lcc_h$group2),  levels = c("HC", "LRRK2_PD"))

ppmi_h$sex2 <- factor(as.character(ppmi_h$sex2), levels = c("female", "male"))
lcc_h$sex2  <- factor(as.character(lcc_h$sex2),  levels = c("female", "male"))

max_missing_train <- 0.30
max_missing_ext   <- 0.50

keep_genes <- common_genes[
  colMeans(is.na(ppmi_h[, common_genes, drop = FALSE])) <= max_missing_train &
    colMeans(is.na(lcc_h[, common_genes, drop = FALSE])) <= max_missing_ext
]

cat("Proteins after missingness filter:", length(keep_genes), "\n")

ppmi_h <- ppmi_h[, c("group2", "age2", "sex2", keep_genes)]
lcc_h  <- lcc_h[,  c("group2", "age2", "sex2", keep_genes)]

cat("Final class counts after harmonization:\n")
print(table(ppmi_h$group2, useNA = "ifany"))
print(table(lcc_h$group2, useNA = "ifany"))
print(table(ppmi_h$sex2, useNA = "ifany"))
print(table(lcc_h$sex2, useNA = "ifany"))

write.csv(
  data.frame(feature = keep_genes),
  file.path(outdir, "shared_proteins_used.csv"),
  row.names = FALSE
)

# =========================
# SHARED INTERNAL CV FOLDS
# =========================
set.seed(23)
nfolds_inner <- 10
fold_id <- sample(rep(1:nfolds_inner, length.out = nrow(ppmi_h)))

# =========================================================
# MODEL 1: ELASTIC NET AGE + SEX ONLY
# =========================================================
oof_pred_cov <- rep(NA_real_, nrow(ppmi_h))

for (k in 1:nfolds_inner) {
  train_idx <- which(fold_id != k)
  val_idx   <- which(fold_id == k)
  
  ppmi_train_fold <- ppmi_h[train_idx, c("group2", "age2", "sex2"), drop = FALSE]
  ppmi_val_fold   <- ppmi_h[val_idx,   c("group2", "age2", "sex2"), drop = FALSE]
  
  imp_fold <- median_impute_by_train(ppmi_train_fold, ppmi_val_fold, "age2")
  train_imp_fold <- imp_fold$train
  val_imp_fold   <- imp_fold$test
  
  sc_fold <- scale_by_train(train_imp_fold, val_imp_fold, "age2")
  train_sc_fold <- sc_fold$train
  val_sc_fold   <- sc_fold$test
  
  x_train_fold <- model.matrix(group2 ~ age2 + sex2, data = train_sc_fold)[, -1, drop = FALSE]
  y_train_fold <- train_sc_fold$group2
  x_val_fold   <- model.matrix(group2 ~ age2 + sex2, data = val_sc_fold)[, -1, drop = FALSE]
  
  cvfit_fold_cov <- cv.glmnet(
    x = x_train_fold,
    y = y_train_fold,
    family = "binomial",
    alpha = 0.5,
    type.measure = "auc",
    nfolds = 10
  )
  
  oof_pred_cov[val_idx] <- as.numeric(
    predict(
      cvfit_fold_cov,
      newx = x_val_fold,
      s = "lambda.min",
      type = "response"
    )
  )
}

res_cov_internal <- calc_auc_ci(ppmi_h$group2, oof_pred_cov)
write.csv(
  res_cov_internal,
  file.path(outdir, "elastic_net_age_sex_ppmi_internal_cv_auc.csv"),
  row.names = FALSE
)

plot_roc_save(
  ppmi_h$group2, oof_pred_cov,
  file.path(outdir, "elastic_net_age_sex_ppmi_internal_cv_ROC.pdf"),
  "Elastic net age + sex | PPMI internal CV"
)

ppmi_cov <- ppmi_h[, c("group2", "age2", "sex2"), drop = FALSE]
lcc_cov  <- lcc_h[,  c("group2", "age2", "sex2"), drop = FALSE]

imp_cov <- median_impute_by_train(ppmi_cov, lcc_cov, "age2")
ppmi_cov_imp <- imp_cov$train
lcc_cov_imp  <- imp_cov$test

sc_cov <- scale_by_train(ppmi_cov_imp, lcc_cov_imp, "age2")
ppmi_cov_sc <- sc_cov$train
lcc_cov_sc  <- sc_cov$test

x_ppmi_cov <- model.matrix(group2 ~ age2 + sex2, data = ppmi_cov_sc)[, -1, drop = FALSE]
y_ppmi_cov <- ppmi_cov_sc$group2

cvfit_cov <- cv.glmnet(
  x = x_ppmi_cov,
  y = y_ppmi_cov,
  family = "binomial",
  alpha = 0.5,
  type.measure = "auc",
  nfolds = 10
)

pred_lcc_cov <- as.numeric(
  predict(
    cvfit_cov,
    newx = model.matrix(group2 ~ age2 + sex2, data = lcc_cov_sc)[, -1, drop = FALSE],
    s = "lambda.min",
    type = "response"
  )
)

res_cov <- calc_auc_ci(lcc_cov_sc$group2, pred_lcc_cov)
write.csv(res_cov, file.path(outdir, "age_sex_only_lcc_auc.csv"), row.names = FALSE)

plot_roc_save(
  lcc_cov_sc$group2, pred_lcc_cov,
  file.path(outdir, "age_sex_only_lcc_ROC.pdf"),
  "Age + sex only | LCC external validation"
)

# =========================================================
# MODEL 2: ELASTIC NET PROTEINS ONLY
# =========================================================
oof_pred <- rep(NA_real_, nrow(ppmi_h))

for (k in 1:nfolds_inner) {
  train_idx <- which(fold_id != k)
  val_idx   <- which(fold_id == k)
  
  ppmi_train_fold <- ppmi_h[train_idx, , drop = FALSE]
  ppmi_val_fold   <- ppmi_h[val_idx, , drop = FALSE]
  
  imp_fold <- median_impute_by_train(ppmi_train_fold, ppmi_val_fold, keep_genes)
  train_imp_fold <- imp_fold$train
  val_imp_fold   <- imp_fold$test
  
  sc_fold <- scale_by_train(train_imp_fold, val_imp_fold, keep_genes)
  train_sc_fold <- sc_fold$train
  val_sc_fold   <- sc_fold$test
  
  x_train_fold <- as.matrix(train_sc_fold[, keep_genes, drop = FALSE])
  y_train_fold <- train_sc_fold$group2
  x_val_fold   <- as.matrix(val_sc_fold[, keep_genes, drop = FALSE])
  
  cvfit_fold <- cv.glmnet(
    x = x_train_fold,
    y = y_train_fold,
    family = "binomial",
    alpha = 0.5,
    type.measure = "auc",
    nfolds = 10
  )
  
  oof_pred[val_idx] <- as.numeric(
    predict(
      cvfit_fold,
      newx = x_val_fold,
      s = "lambda.min",
      type = "response"
    )
  )
}

res_enet_internal <- calc_auc_ci(ppmi_h$group2, oof_pred)
write.csv(
  res_enet_internal,
  file.path(outdir, "elastic_net_ppmi_internal_cv_auc.csv"),
  row.names = FALSE
)

plot_roc_save(
  ppmi_h$group2, oof_pred,
  file.path(outdir, "elastic_net_ppmi_internal_cv_ROC.pdf"),
  "Elastic net proteins only | PPMI internal CV"
)

imp <- median_impute_by_train(ppmi_h, lcc_h, keep_genes)
ppmi_imp <- imp$train
lcc_imp  <- imp$test

sc <- scale_by_train(ppmi_imp, lcc_imp, keep_genes)
ppmi_sc <- sc$train
lcc_sc  <- sc$test

x_ppmi <- as.matrix(ppmi_sc[, keep_genes, drop = FALSE])
y_ppmi <- ppmi_sc$group2

cvfit <- cv.glmnet(
  x = x_ppmi,
  y = y_ppmi,
  family = "binomial",
  alpha = 0.5,
  type.measure = "auc",
  nfolds = 10
)

coef_df <- data.frame(
  feature = rownames(as.matrix(coef(cvfit, s = "lambda.min"))),
  coef = as.numeric(coef(cvfit, s = "lambda.min"))
)
coef_df <- subset(coef_df, coef != 0)
write.csv(
  coef_df,
  file.path(outdir, "elastic_net_selected_features.csv"),
  row.names = FALSE
)

pred_lcc_enet <- as.numeric(
  predict(
    cvfit,
    newx = as.matrix(lcc_sc[, keep_genes, drop = FALSE]),
    s = "lambda.min",
    type = "response"
  )
)

res_enet <- calc_auc_ci(lcc_sc$group2, pred_lcc_enet)
write.csv(
  res_enet,
  file.path(outdir, "elastic_net_lcc_auc.csv"),
  row.names = FALSE
)

plot_roc_save(
  lcc_sc$group2, pred_lcc_enet,
  file.path(outdir, "elastic_net_lcc_ROC.pdf"),
  "Elastic net proteins only | LCC external validation"
)

# =========================================================
# MODEL 3: ELASTIC NET PROTEINS + AGE + SEX
# =========================================================
oof_pred_full <- rep(NA_real_, nrow(ppmi_h))

for (k in 1:nfolds_inner) {
  train_idx <- which(fold_id != k)
  val_idx   <- which(fold_id == k)
  
  ppmi_train_fold <- ppmi_h[train_idx, c("group2", "age2", "sex2", keep_genes), drop = FALSE]
  ppmi_val_fold   <- ppmi_h[val_idx,   c("group2", "age2", "sex2", keep_genes), drop = FALSE]
  
  imp_fold <- median_impute_by_train(ppmi_train_fold, ppmi_val_fold, c("age2", keep_genes))
  train_imp_fold <- imp_fold$train
  val_imp_fold   <- imp_fold$test
  
  sc_fold <- scale_by_train(train_imp_fold, val_imp_fold, c("age2", keep_genes))
  train_sc_fold <- sc_fold$train
  val_sc_fold   <- sc_fold$test
  
  x_train_fold <- model.matrix(
    group2 ~ age2 + sex2 + .,
    data = train_sc_fold[, c("group2", "age2", "sex2", keep_genes), drop = FALSE]
  )[, -1, drop = FALSE]
  
  y_train_fold <- train_sc_fold$group2
  
  x_val_fold <- model.matrix(
    group2 ~ age2 + sex2 + .,
    data = val_sc_fold[, c("group2", "age2", "sex2", keep_genes), drop = FALSE]
  )[, -1, drop = FALSE]
  
  cvfit_fold_full <- cv.glmnet(
    x = x_train_fold,
    y = y_train_fold,
    family = "binomial",
    alpha = 0.5,
    type.measure = "auc",
    nfolds = 10
  )
  
  oof_pred_full[val_idx] <- as.numeric(
    predict(
      cvfit_fold_full,
      newx = x_val_fold,
      s = "lambda.min",
      type = "response"
    )
  )
}

res_enet_full_internal <- calc_auc_ci(ppmi_h$group2, oof_pred_full)
write.csv(
  res_enet_full_internal,
  file.path(outdir, "elastic_net_age_sex_proteins_ppmi_internal_cv_auc.csv"),
  row.names = FALSE
)

plot_roc_save(
  ppmi_h$group2, oof_pred_full,
  file.path(outdir, "elastic_net_age_sex_proteins_ppmi_internal_cv_ROC.pdf"),
  "Elastic net age + sex + proteins | PPMI internal CV"
)

ppmi_full_enet <- ppmi_h[, c("group2", "age2", "sex2", keep_genes), drop = FALSE]
lcc_full_enet  <- lcc_h[,  c("group2", "age2", "sex2", keep_genes), drop = FALSE]

imp_full <- median_impute_by_train(ppmi_full_enet, lcc_full_enet, c("age2", keep_genes))
ppmi_full_imp <- imp_full$train
lcc_full_imp  <- imp_full$test

sc_full <- scale_by_train(ppmi_full_imp, lcc_full_imp, c("age2", keep_genes))
ppmi_full_sc <- sc_full$train
lcc_full_sc  <- sc_full$test

x_ppmi_full <- model.matrix(
  group2 ~ age2 + sex2 + .,
  data = ppmi_full_sc
)[, -1, drop = FALSE]

y_ppmi_full <- ppmi_full_sc$group2

cvfit_full_enet <- cv.glmnet(
  x = x_ppmi_full,
  y = y_ppmi_full,
  family = "binomial",
  alpha = 0.5,
  type.measure = "auc",
  nfolds = 10
)

coef_df_full <- data.frame(
  feature = rownames(as.matrix(coef(cvfit_full_enet, s = "lambda.min"))),
  coef = as.numeric(coef(cvfit_full_enet, s = "lambda.min"))
)
coef_df_full <- subset(coef_df_full, coef != 0)
write.csv(
  coef_df_full,
  file.path(outdir, "elastic_net_age_sex_proteins_selected_features.csv"),
  row.names = FALSE
)

x_lcc_full <- model.matrix(
  group2 ~ age2 + sex2 + .,
  data = lcc_full_sc
)[, -1, drop = FALSE]

pred_lcc_enet_full <- as.numeric(
  predict(
    cvfit_full_enet,
    newx = x_lcc_full,
    s = "lambda.min",
    type = "response"
  )
)

res_enet_full <- calc_auc_ci(lcc_full_sc$group2, pred_lcc_enet_full)
write.csv(
  res_enet_full,
  file.path(outdir, "elastic_net_age_sex_proteins_lcc_auc.csv"),
  row.names = FALSE
)

plot_roc_save(
  lcc_full_sc$group2, pred_lcc_enet_full,
  file.path(outdir, "elastic_net_age_sex_proteins_lcc_ROC.pdf"),
  "Elastic net age + sex + proteins | LCC external validation"
)

# =========================
# 12. SUMMARY FILE
# =========================
summary_df <- rbind(
  cbind(model = "Age + sex only", res_cov),
  cbind(model = "Elastic net proteins only", res_enet),
  cbind(model = "Elastic net age + sex + proteins", res_enet_full)
)

write.csv(
  summary_df,
  file.path(outdir, "model_summary_lcc_external_validation.csv"),
  row.names = FALSE
)

# =========================
# MAIN FIGURE:
# 3 external ROC curves
# =========================
roc_cov <- roc(
  response = lcc_cov_sc$group2,
  predictor = pred_lcc_cov,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

roc_enet <- roc(
  response = lcc_sc$group2,
  predictor = pred_lcc_enet,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

roc_enet_full <- roc(
  response = lcc_full_sc$group2,
  predictor = pred_lcc_enet_full,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

pdf(file.path(outdir, "combined_ROC_3models_external.pdf"), width = 6, height = 6)

plot.roc(
  roc_cov,
  legacy.axes = TRUE,
  col = "#1b9e77",
  lwd = 2,
  print.auc = FALSE,
  main = "LCC external validation",
  xlab = "1 - Specificity",
  ylab = "Sensitivity"
)

plot.roc(
  roc_enet,
  legacy.axes = TRUE,
  col = "#d95f02",
  lwd = 2,
  print.auc = FALSE,
  add = TRUE
)

plot.roc(
  roc_enet_full,
  legacy.axes = TRUE,
  col = "#7570b3",
  lwd = 2,
  print.auc = FALSE,
  add = TRUE
)

abline(a = 0, b = 1, lty = 2, col = "grey")

legend(
  "bottomright",
  legend = c(
    paste0("Age + sex only (AUC=", sprintf("%.3f", as.numeric(auc(roc_cov))), ")"),
    paste0("Proteins only (AUC=", sprintf("%.3f", as.numeric(auc(roc_enet))), ")"),
    paste0("Proteins + age + sex (AUC=", sprintf("%.3f", as.numeric(auc(roc_enet_full))), ")")
  ),
  col = c("#1b9e77", "#d95f02", "#7570b3"),
  lwd = 2,
  bty = "n"
)

dev.off()

# =========================
# SUPPLEMENTARY FIGURE:
# 3 internal CV ROC curves
# =========================
roc_cov_internal <- roc(
  response = ppmi_h$group2,
  predictor = oof_pred_cov,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

roc_enet_internal <- roc(
  response = ppmi_h$group2,
  predictor = oof_pred,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

roc_enet_full_internal <- roc(
  response = ppmi_h$group2,
  predictor = oof_pred_full,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

pdf(file.path(outdir, "combined_ROC_3models_internalCV.pdf"), width = 6, height = 6)

plot.roc(
  roc_cov_internal,
  legacy.axes = TRUE,
  col = "#1b9e77",
  lwd = 2,
  print.auc = FALSE,
  main = "PPMI internal CV",
  xlab = "1 - Specificity",
  ylab = "Sensitivity"
)

plot.roc(
  roc_enet_internal,
  legacy.axes = TRUE,
  col = "#d95f02",
  lwd = 2,
  print.auc = FALSE,
  add = TRUE
)

plot.roc(
  roc_enet_full_internal,
  legacy.axes = TRUE,
  col = "#7570b3",
  lwd = 2,
  print.auc = FALSE,
  add = TRUE
)

abline(a = 0, b = 1, lty = 2, col = "grey")

legend(
  "bottomright",
  legend = c(
    paste0("Age + sex only (AUC=", sprintf("%.3f", as.numeric(auc(roc_cov_internal))), ")"),
    paste0("Proteins only (AUC=", sprintf("%.3f", as.numeric(auc(roc_enet_internal))), ")"),
    paste0("Proteins + age + sex (AUC=", sprintf("%.3f", as.numeric(auc(roc_enet_full_internal))), ")")
  ),
  col = c("#1b9e77", "#d95f02", "#7570b3"),
  lwd = 2,
  bty = "n"
)

dev.off()

# =========================
# EXTRA PLOT:
# proteins-only external ROC
# =========================
roc_enet_clean <- roc(
  response = lcc_sc$group2,
  predictor = pred_lcc_enet,
  levels = c("HC", "LRRK2_PD"),
  direction = "<"
)

pdf(file.path(outdir, "elastic_net_lcc_ROC_clean.pdf"), width = 5, height = 5)
plot.roc(
  roc_enet_clean,
  legacy.axes = TRUE,
  print.auc = TRUE,
  col = "black",
  lwd = 2,
  main = "Elastic net: LRRK2 PD vs HC\nExternal validation in LCC"
)
dev.off()

# =========================
# PR FIGURE + AP SUMMARY
# best external model
# default set to proteins-only elastic net
# If proteins+age+sex is better, replace lcc_sc/pred_lcc_enet with lcc_full_sc/pred_lcc_enet_full
# =========================
pr_dat <- data.frame(
  truth = factor(as.character(lcc_sc$group2), levels = c("HC", "LRRK2_PD")),
  .pred_LRRK2_PD = pred_lcc_enet
)

ap_res <- average_precision(
  pr_dat,
  truth = truth,
  .pred_LRRK2_PD,
  event_level = "second"
)
write.csv(ap_res, file.path(outdir, "elastic_net_lcc_average_precision.csv"), row.names = FALSE)

pr_curve_df <- pr_curve(
  pr_dat,
  truth = truth,
  .pred_LRRK2_PD,
  event_level = "second"
)
write.csv(pr_curve_df, file.path(outdir, "elastic_net_lcc_PR_points.csv"), row.names = FALSE)

pdf(file.path(outdir, "elastic_net_lcc_PR_curve.pdf"), width = 5, height = 5)
print(
  autoplot(pr_curve_df) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste0(
        "Elastic net PR curve: LRRK2 PD vs HC\nExternal validation in LCC (AP = ",
        round(ap_res$.estimate, 3), ")"
      ),
      x = "Recall",
      y = "Precision"
    )
)
dev.off()

# =========================
# ELASTIC NET FEATURE PLOT
# =========================
coef_plot_df <- coef_df
coef_plot_df <- coef_plot_df[coef_plot_df$feature != "(Intercept)", , drop = FALSE]
coef_plot_df$abs_coef <- abs(coef_plot_df$coef)
coef_plot_df <- coef_plot_df[order(coef_plot_df$abs_coef, decreasing = TRUE), , drop = FALSE]

top_n <- min(20, nrow(coef_plot_df))
coef_top <- coef_plot_df[1:top_n, , drop = FALSE]
coef_top$feature <- factor(coef_top$feature, levels = rev(coef_top$feature))

pdf(file.path(outdir, "elastic_net_top_features.pdf"), width = 6, height = 7)
print(
  ggplot(coef_top, aes(x = feature, y = coef)) +
    geom_col(fill = "#f5535d") +
    coord_flip() +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    ) +
    labs(
      title = "Top elastic-net features",
      x = "",
      y = "Coefficient"
    )
)
dev.off()

# =========================
# MODEL COMPARISON PLOT
# =========================
summary_df$AUC <- as.numeric(summary_df$AUC)
summary_df$CI_lower <- as.numeric(summary_df$CI_lower)
summary_df$CI_upper <- as.numeric(summary_df$CI_upper)

summary_df$model <- factor(
  summary_df$model,
  levels = c(
    "Age + sex only",
    "Elastic net proteins only",
    "Elastic net age + sex + proteins"
  )
)

pdf(file.path(outdir, "model_comparison_AUC_barplot.pdf"), width = 7, height = 5)
print(
  ggplot(summary_df, aes(x = model, y = AUC)) +
    geom_col(fill = "#f5535d") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    theme_bw(base_size = 13) +
    ylim(0, 1) +
    labs(
      title = "External validation in LCC",
      x = "",
      y = "ROC AUC"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
)
dev.off()

cat("Done.\n")
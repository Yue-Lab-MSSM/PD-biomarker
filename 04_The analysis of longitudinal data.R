library(lme4)
library(tidyverse)
library(afex)
library(dplyr)
# Example data preparation (replace with your actual data loading)
df <- read.csv("project_177_csf without ST_with years after onset.csv")
df <- df %>% filter(SexM1.F0 == 1)
#df <- df %>% filter(SexM1.F0 == 1 & EVENT_ID <= 5)


df$time <- df$onset_EVENT_ID
# Check for unmapped EVENT_ID values
if (any(is.na(df$time))) {
  warning("Some EVENT_ID values could not be mapped to time. Check your data.")
}

# View the first few rows to confirm
head(df[, c("PATNO", "onset_EVENT_ID", "time")])


####Count Unique Time Points per Participant >2
# Step 1: Count the number of unique time points per PATNO
visit_counts <- df %>%
  group_by(PATNO) %>%
  summarise(n_visits = n_distinct(time)) %>%
  ungroup()

# Step 2: Identify PATNOs with more than 2 visits
valid_patnos <- visit_counts %>%
  filter(n_visits > 0) %>%
  pull(PATNO)

# Step 3: Filter the original dataframe to keep only these participants
df <- df %>%
  filter(PATNO %in% valid_patnos)

# Step 4: Verify the result
visit_counts_filtered <- df %>%
  group_by(PATNO) %>%
  summarise(n_visits = n_distinct(time)) %>%
  ungroup()
print(visit_counts_filtered)



# Define group labels (example)
df <- df %>%
  mutate(group_label = case_when(
    Group == "HC" ~ "HC",
    Group == "PD_LRRK2" ~ "PD_LRRK2",
    TRUE ~ "Other"
  ))




# Filter for relevant groups
long_df <- df %>% filter(group_label %in% c("HC", "PD_LRRK2"))

# Identify protein columns (adjust the pattern based on your column names)
protein_cols <- grep("^P177_CSF_", names(long_df), value = TRUE)

# Apply log2 transformation to all protein columns
long_df <- long_df %>%
  mutate(across(all_of(protein_cols), ~ log2(. + 1), .names = "{.col}"))



######################################
#long_df <- long_df[!(long_df$PATNO %in% c(" ")), ]  #One participant is a outlier, as it shows >2SD of intercept for many proteins, because of the PPMI policy, we can not show the participant ID
#######################################


# Verify the transformation by checking a few rows
head(long_df[, c("PATNO", "time", protein_cols[1:3])])

# Initialize a list to store results
results <- list()


##Principal Component Analysis (PCA)
# 2. Subset your data frame to just those columns
df_prot <- long_df[, protein_cols]


# 2. Impute via missMDA
library(missMDA)
imp      <- imputePCA(df_prot, ncp = 3, method = "EM")
prot_imp <- imp$completeObs

# 3a. (Preferred) Run PCA directly: rows = samples, cols = proteins
pc <- prcomp(prot_imp, center = TRUE, scale. = TRUE)

library(ggplot2)

# Build a data.frame of PCs and metadata
pc_df <- data.frame(PC1 = pc$x[,1],
                    PC2 = pc$x[,2],
                    Group = long_df$Group,
                    PATNO = long_df$sudo.PATNO)

ggplot(pc_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 1) +
  geom_text(aes(label = PATNO), 
            check_overlap = TRUE, 
            size = 3, 
            alpha = 0.6) +
  labs(title = "PCA of CSF Proteome",
       x = "PC1", y = "PC2") +
  theme_minimal()



# 2. Inspect age distribution
summary(long_df$envolved.age)

table(long_df$Group)

long_df %>%
  distinct(PATNO, group_label) %>%
  count(group_label)

# 3. Mean‐center envolved.age
long_df$envolved.age_c <- long_df$envolved.age - mean(long_df$envolved.age, na.rm = TRUE)


library(lme4)
library(lmerTest)
library(dplyr)

# Containers
all_AICs      <- list()
all_summaries <- list()

# 1. Model selection (ML) & record AICs
for (protein in protein_cols) {
  
  # Define candidate formulas (as before)
  formula_reduced_slope <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time + group_label:time + (1 + time | PATNO)")
  )
  formula_reduced_non_random <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time + group_label:time + (1 | PATNO)")
  )
  formula_full_slope <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time + group_label:time + envolved.age_c + (1 + time | PATNO)")
  )
  formula_full_non_random <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time + group_label:time + envolved.age_c + (1 | PATNO)")
  )
  
  # 1a. Fit A2 (slope) under ML
  fit_reduced_slope_ml <- tryCatch(
    lmer(formula_reduced_slope, data = long_df, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) e
  )
  if (inherits(fit_reduced_slope_ml, "error") || isSingular(fit_reduced_slope_ml, tol = 1e-4)) {
    fit_A_final_ml <- tryCatch(
      lmer(formula_reduced_non_random, data = long_df, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) e
    )
    reason_A <- if (inherits(fit_reduced_slope_ml, "error")) "error_reduced_slope_REML"
    else paste0("reduced_slope_non_random_REML")
  } else {
    fit_A_final_ml <- fit_reduced_slope_ml
    reason_A <- "reduced_slope_REML"
  }
  if (!inherits(fit_A_final_ml, "error") && isSingular(fit_A_final_ml, tol = 1e-4)) {
    reason_A <- paste0(reason_A, "_reduced_non_random_REML")
  }
  
  # 1b. Fit B2 (full) under ML
  fit_full_slope_ml <- tryCatch(
    lmer(formula_full_slope, data = long_df, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) e
  )
  if (inherits(fit_full_slope_ml, "error") || isSingular(fit_full_slope_ml, tol = 1e-4)) {
    fit_full_final_ml <- tryCatch(
      lmer(formula_full_non_random, data = long_df, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) e
    )
    reason_full <- if (inherits(fit_full_slope_ml, "error")) "error_full_slope_REML"
    else paste0("full_slope_non_random_REML")
  } else {
    fit_full_final_ml <- fit_full_slope_ml
    reason_full <- "full_slope_REML"
  }
  if (!inherits(fit_full_final_ml, "error") && isSingular(fit_full_final_ml, tol = 1e-4)) {
    reason_full <- paste0(reason_full, "_full_non_random_REML")
  }
  
  # 1c. If both fail, skip
  if (inherits(fit_A_final_ml, "error") && inherits(fit_full_final_ml, "error")) {
    warning("Protein '", protein, "' could not be fitted under REML.")
    next
  }
  
  # 1d. Record AICs (for comparison under ML)
  AIC_A_ml    <- if (!inherits(fit_A_final_ml, "error"))    AIC(fit_A_final_ml)    else Inf
  AIC_full_ml <- if (!inherits(fit_full_final_ml, "error")) AIC(fit_full_final_ml) else Inf
  
  if (AIC_full_ml < AIC_A_ml) {
    chosen_model_ml       <- fit_full_final_ml
    chosen_fixed_spec_ml  <- "with_envolved.age_c"
    chosen_random_ml      <- reason_full
  } else {
    chosen_model_ml       <- fit_A_final_ml
    chosen_fixed_spec_ml  <- "without_envolved.age_c"
    chosen_random_ml      <- reason_A
  }
  
  all_AICs[[protein]] <- data.frame(
    protein            = protein,
    AIC_without_ENROLL = AIC_A_ml,
    AIC_with_ENROLL    = AIC_full_ml,
    chosen_fixed_ml    = chosen_fixed_spec_ml,
    chosen_random_ml   = chosen_random_ml,
    stringsAsFactors   = FALSE
  )
  
  
  # 2. Refit the chosen model under REML for final estimates
  #    We extract the formula from the ML‐fit object and simply set REML = TRUE
  best_formula <- formula(chosen_model_ml)
  fit_reml_final <- update(chosen_model_ml, REML = TRUE)
  
  # 2a. Store the summary of the REML‐fit
  summ_reml <- summary(fit_reml_final)
  coef_mat  <- summ_reml$coefficients
  
  capture.output(
    summ_reml,
    file = paste0(
      protein,
      "_summary_of_fix_effect_random.txt"))
  
  
  # 2b. Extract relevant fields (adjust as needed)
  summary_row <- data.frame(
    protein                  = protein,
    chosen_fixed_REML        = chosen_fixed_spec_ml,
    chosen_random_REML       = chosen_random_ml,
    est_group_labelHC_REML   = if("(Intercept)"   %in% rownames(coef_mat)) coef_mat["(Intercept)",   "Estimate"] else NA,
    pval_group_labelHC_REML  = if("(Intercept)"   %in% rownames(coef_mat)) coef_mat["(Intercept)",   "Pr(>|t|)"]  else NA,
    est_time_REML            = if("time"             %in% rownames(coef_mat)) coef_mat["time",             "Estimate"] else NA,
    pval_time_REML           = if("time"             %in% rownames(coef_mat)) coef_mat["time",             "Pr(>|t|)"]  else NA,
    est_group_PD_LRRK2_REML  = if("group_labelPD_LRRK2" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2", "Estimate"] else NA,
    pval_group_PD_LRRK2_REML = if("group_labelPD_LRRK2" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2", "Pr(>|t|)"]  else NA,
    est_group_labelPD_LRRK2_time_REML     = if("group_labelPD_LRRK2:time" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2:time", "Estimate"] else NA,
    pval_group_labelPD_LRRK2_time_REML    = if("group_labelPD_LRRK2:time" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2:time", "Pr(>|t|)"]  else NA,
    est_envolved.age_c_REML    = if("envolved.age_c"     %in% rownames(coef_mat)) coef_mat["envolved.age_c",     "Estimate"] else NA,
    pval_envolved.age_c_REML   = if("envolved.age_c"     %in% rownames(coef_mat)) coef_mat["envolved.age_c",     "Pr(>|t|)"]  else NA,
    var_PATNO_REML           = as.data.frame(summ_reml$varcor)$vcov[1],
    sd_PATNO_REML            = as.data.frame(summ_reml$varcor)$sdcor[1],
    n_obs_REML               = summ_reml$devcomp$dims["n"],
    n_groups_REML            = summ_reml$ngrps[["PATNO"]],
    stringsAsFactors         = FALSE
  )
  all_summaries[[protein]] <- summary_row
}  

# 3. Combine everything at the end and write CSVs
df_AICs      <- bind_rows(all_AICs)
df_summaries <- bind_rows(all_summaries)

write.csv(df_AICs,      "all_proteins_AIC_comparisons_ML.csv",    row.names = FALSE)
write.csv(df_summaries, "all_proteins_model_summaries_REML.csv", row.names = FALSE)


##Volcano Plot of Effect Sizes
library(EnhancedVolcano)

library(dplyr)

df_summaries <- bind_rows(all_summaries)
volcano_df <- df_summaries %>%
  # Keep only proteins for which `est_group_PD_LRRK2_REML` is not NA:
  filter(!is.na(est_group_PD_LRRK2_REML) & !is.na(pval_group_PD_LRRK2_REML)) %>%
  select(
    Protein  = protein,
    log2FC   = est_group_PD_LRRK2_REML,
    pvalue   = pval_group_PD_LRRK2_REML
  ) %>%
  # (Optional) If your protein column has a prefix like "P177_CSF_", remove it:
  mutate(
    Protein = sub("^([^_]*_){2}", "", Protein)  # keep text after second underscore
  ) %>%
  # Precompute –log10(pvalue)
  mutate(
    negLog10P = -log10(pvalue)
  )# handy to precompute


pdf("The volcano plot of BL_FC and p value for PB.pdf",7.5,4.5)
EnhancedVolcano(volcano_df,
                #selectLab = " ",
                lab           = volcano_df$Protein,               # name of the column with labels
                x             = 'log2FC',       # the log₂FC column
                y             = 'pvalue',    # the p-value column
                xlab          = bquote(~Log[2]~ "FC"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.05,                     # raw p-value cutoff
                FCcutoff      = 0.1,
                xlim = c(-0.5,0.5),
                ylim = c(0,2),
                legendLabels=c('NS','Log2FC','pvalue','Log2FC & pvalue'),
                gridlines.major = F, gridlines.minor = F,
                pointSize = 1,
                legendPosition = "right",
                legendLabSize = 14,
                legendIconSize = 2,
                legendDropLevels = F,
                #colAlpha = 0.5,
                drawConnectors = T, min.segment.length = 0.0001, 
                selectLab = c("CTSB", "CTSD", "NPC2", "PSAP", "ICOSLG", "B2M", "SERPING1", "CD14",
                              "CLU","COL6A1","LAMP2","TGFBI","EFEMP1", "FBLN5","APOE", "FN1", "CSF1R", "CHI3L1"),
                title = NULL,subtitle = "LRRK2_PD vs. HC male",
                max.overlaps = 10
)

dev.off()


volcano_df1 <- df_summaries %>%
  filter(!is.na(est_group_labelPD_LRRK2_time_REML) & !is.na(pval_group_labelPD_LRRK2_time_REML)) %>%
  select(
    Protein  = protein,
    log2FC   = est_group_labelPD_LRRK2_time_REML,
    pvalue   = pval_group_labelPD_LRRK2_time_REML
  ) %>%
  # (Optional) If your protein column has a prefix like "P177_CSF_", remove it:
  mutate(
    Protein = sub("^([^_]*_){2}", "", Protein)  # keep text after second underscore
  ) %>%
  # Precompute –log10(pvalue)
  mutate(
    negLog10P = -log10(pvalue)
  )# handy to precompute



pdf("The volcano plot of disease slope_change and p value for PB.pdf",7.5,4.5)
EnhancedVolcano(volcano_df1,
                lab           = volcano_df1$Protein,               # name of the column with labels
                x             = 'log2FC',       # the log₂FC column
                y             = 'pvalue',    # the p-value column
                xlab          = bquote(~"Disease_Slope_Change"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.05,                     # raw p-value cutoff
                FCcutoff      = 0.01,
                legendLabels=c('NS','Log2FC','pvalue','Log2FC & pvalue'),
                gridlines.major = F, gridlines.minor = F,
                xlim = c(-0.1,0.1),
                ylim = c(0,3),
                pointSize = 1,
                legendPosition = "right",
                legendLabSize = 14,
                legendIconSize = 2,
                legendDropLevels = F,
                colAlpha = 0.5,
                drawConnectors = T, min.segment.length = 0.0001, 
                selectLab = c("CTSB", "CTSD", "NPC2", "PSAP", "ICOSLG", "B2M", "SERPING1", "CD14",
                              "CLU","COL6A1","LAMP2","TGFBI","EFEMP1", "FBLN5","APOE", "FN1", "CSF1R", "CHI3L1"),
                title = NULL,subtitle = "LRRK2_PD vs. HC male",
                max.overlaps = 10
)


dev.off()


volcano_df2 <- df_summaries %>%
  filter(!is.na(est_time_REML) & !is.na(pval_time_REML)) %>%
  select(
    Protein  = protein,
    log2FC   = est_time_REML,
    pvalue   = pval_time_REML
  ) %>%
  # (Optional) If your protein column has a prefix like "P177_CSF_", remove it:
  mutate(
    Protein = sub("^([^_]*_){2}", "", Protein)  # keep text after second underscore
  ) %>%
  # Precompute –log10(pvalue)
  mutate(
    negLog10P = -log10(pvalue)
  )# handy to precompute


pdf("The volcano plot of HC slope_change and p value for PB.pdf",7.5,4.5)
EnhancedVolcano(volcano_df2,
                #selectLab = " ",
                lab           = volcano_df2$Protein,               # name of the column with labels
                x             = 'log2FC',       # the log₂FC column
                y             = 'pvalue',    # the p-value column
                xlab          = bquote(~"HC Slope_Change"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.05,                     # raw p-value cutoff
                FCcutoff      = 0.01,
                xlim = c(-0.1,0.1),
                ylim = c(0,2),
                legendLabels=c('NS','Log2FC','pvalue','Log2FC & pvalue'),
                gridlines.major = F, gridlines.minor = F,
                pointSize = 1,
                legendPosition = "right",
                legendLabSize = 14,
                legendIconSize = 2,
                legendDropLevels = F,
                #colAlpha = 0.5,
                drawConnectors = T, min.segment.length = 0.0001, 
                selectLab = c(" "),
                title = NULL,subtitle = "LRRK2_PD vs. HC male",
                max.overlaps = 10
)

dev.off()





library(ggplot2)
library(lattice)
library(ggeffects)
library(effects)
library(broom)    # for tidy()
library(dplyr)
library(ggplot2)
library(rlang)
library(broom.mixed)
# 1. Ensure Group is a factor

long_df$group_label <- factor(long_df$group_label)


protein_name <- c("CTSD","LAMP2","GM2A","CLU","CD14","B2M","EFEMP1","APOE","FN1","NPC2","PSAP","CTSB",
                  "ICOSLG","SERPING1","CSF1R","CHI3L1","TGFBI","FBLN5")
proteins  <- paste0("P177_CSF_", protein_name)

# 2. Fit by REML
for (response_var in proteins)  {
  
  fml <- as.formula(
    paste0(
      response_var,
      " ~ 1 + group_label + time + group_label:time +   (1  + time | PATNO)"      # envolved.age_c +       + time
    )
  )
  
  m_mixed <- mixed(
    formula = fml,
    data    = long_df,
    control = lmerControl(optimizer = "bobyqa"),
    method  = "LRT"
  )
  
  capture.output(
    m_mixed,
    file = paste0(response_var, "_fix_effect_significance_test_by_mix.txt")
  )
  
  m_reml <- lmer(
    formula = fml,
    data    = long_df,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  capture.output(
    summary(m_reml),
    file = paste0(response_var, " _lmer_full_summary.txt")
  )
  
  # 3. Build a prediction grid (marginal: fixed effects only)
  pred_df <- expand.grid(
    time        = seq(min(long_df$time), max(long_df$time), length.out = 100),
    group_label = levels(long_df$group_label)
  )
  pred_df$group_label   <- factor(pred_df$group_label, levels = levels(long_df$group_label))
  pred_df$envolved.age_c <- 0
  
  # 4. Predict marginal fit + 95% CI for fixed effects
  pred_df$fit <- predict(m_reml, newdata = pred_df, re.form = NA)
  
  X   <- model.matrix(~ 1 + group_label + time + group_label:time , data = pred_df)
  V   <- vcov(m_reml)                       # Var-cov of fixed effects
  se  <- sqrt(diag(X %*% V %*% t(X)))       # SE for each row's linear predictor
  
  pred_df$lo <- pred_df$fit - 1.96 * se
  pred_df$hi <- pred_df$fit + 1.96 * se
  
  # 5. (Optional) tidy fixed effects for slopes/p-values used in annotations
  tidy_fix <- broom::tidy(m_reml, effects = "fixed")
  slope_hc <- tidy_fix  %>% filter(term == "time") %>% pull(estimate)
  p_hc     <- tidy_fix  %>% filter(term == "time") %>% pull(p.value)
  
  b_int    <- tidy_fix  %>% filter(term == "group_labelPD_LRRK2:time") %>% pull(estimate)
  p_int    <- tidy_fix  %>% filter(term == "group_labelPD_LRRK2:time") %>% pull(p.value)
  
  slope_pd <- slope_hc + b_int
  p_pd     <- p_int
  
  yrange   <- range(long_df[[response_var]], na.rm = TRUE)
  xpos     <- max(long_df$time) * 0.7
  annot_df <- tibble(
    group_label = c("HC", "PD_LRRK2"),
    x           = xpos,
    y           = c(yrange[2]*0.9, yrange[2]*0.8),
    label       = c(
      sprintf("slope=%.3f\np=%.3g",  slope_hc, p_hc),
      sprintf("slope=%.3f\np=%.3g",  slope_pd, p_pd)
    )
  )
  
  # 6. Plot
  pdf(paste0(response_var, "_ObservedData_Fit_connecting.pdf"), 5, 3)
  tryCatch({
    y_rng    <- range(long_df[[response_var]], na.rm = TRUE)
    y_min    <- floor(y_rng[1])
    y_max    <- ceiling(y_rng[2])
    y_breaks <- seq(y_min, y_max, by = 1)
    
    p <- ggplot(long_df, aes(time, !!sym(response_var), colour = group_label)) +
      geom_point(alpha = 0.3, size = 0.2) +
      
      # 95% CI ribbons per group (marginal fits)
      geom_ribbon(
        data = pred_df,
        aes(x = time, ymin = lo, ymax = hi, fill = group_label),
        alpha = 0.18,
        inherit.aes = FALSE
      ) +
      geom_line(
        data = pred_df,
        aes(time, fit, colour = group_label, group = group_label),
        size = 1
      ) +
      
      geom_line(aes(group = PATNO), alpha = 0.25, size = 0.3)+
      
      # optional: add text annotations for slopes/p-values
      # geom_text(
      #   data = annot_df,
      #   aes(x = x, y = y, label = label, color = group_label),
      #   hjust = 0, vjust = -0.5, size = 3.5, show.legend = FALSE
      # ) +
      
      scale_color_manual(values = c(HC = "black", PD_LRRK2 = "red")) +
      scale_fill_manual(values  = c(HC = "black", PD_LRRK2 = "red")) +
      
      scale_y_continuous(breaks = y_breaks, limits = c(y_min, y_max)) +
      labs(
        x     = "Years",
        y     = response_var,
        title = "Observed Data and Estimated Slopes",
        color = "group",
        fill  = "group"
      ) +
      
      # Clean look: no grid/background
      theme_classic() +
      theme(
        axis.text    = element_text(color = "black", size = 18),
        axis.title   = element_text(color = "black", size = 18),
        legend.title = element_text(color = "black", size = 18),
        legend.text  = element_text(color = "black", size = 18)
      )
    
    print(p)
  }, finally = dev.off())
}


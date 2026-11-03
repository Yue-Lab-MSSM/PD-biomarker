library(lme4)
library(tidyverse)
library(afex)
library(dplyr)
# Example data preparation (replace with your actual data loading)
df <- read.csv("project_177_csf.csv")
df <- df %>% filter(SexM1.F0 == 1) # for female, SexM1.F0 == 0

df$time_since_onset <- df$onset_EVENT_ID
# Check for unmapped EVENT_ID values
if (any(is.na(df$time_since_onset))) {
  warning("Some EVENT_ID values could not be mapped to time_since_onset. Check your data.")
}

# View the first few rows to confirm
head(df[, c("PATNO", "onset_EVENT_ID", "time_since_onset")])


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


long_df <- long_df %>% 
  
  # Apply log2 transformation to all protein columns  
  mutate(across(all_of(protein_cols), ~ log2(. + 1), .names = "{.col}")) %>%
  ## 1 ─ centre covariates that go into interactions
  mutate(
    envolved_age_c = envolved_age - mean(envolved_age, na.rm = TRUE)
  ) %>% 
  
  ## 2 ─ GRAND-MEAN centre every protein (suffix “_gc”)
  mutate(across(all_of(protein_cols),
                ~ .x - mean(.x, na.rm = TRUE),
                .names = "{.col}_gc")) %>% 
  
  ## 3 ─ PERSON means and within-person residuals
  group_by(PATNO) %>% 
  mutate(across(ends_with("_gc"),
                list(between =  ~ mean(.x, na.rm = TRUE),      # person mean
                     within  =  ~ .x - mean(.x, na.rm = TRUE)  # deviation
                ),
                .names = "{.fn}_{.col}")) %>% 
  ungroup() %>% 
  
  ## 4 ─ z-scale the four new terms per protein
  mutate(across(
    ends_with("_gc"),            # picks up P177_CSF_…_gc, between_…_gc, within_…_gc
    ~ as.numeric(scale(.x)),     # z-score each one
    .names = "{.col}_z"
  ))


head(long_df[, c("PATNO", "time_since_onset", protein_cols[1:3])])

# Initialize a list to store results
results <- list()


##Principal Component Analysis (PCA)
# 1. Subset the data frame to just those columns
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


#remove the outliers
long_df <- long_df[!(long_df$sudo.PATNO %in% c("3004_V10", "3636_BL")), ] #outliers by pca
long_df <- long_df[!(long_df$PATNO %in% c("40585")), ]  #40585 is a outlier, as it shows > 2SD of intercept

# Inspect age distribution
summary(long_df$envolved_age)

table(long_df$Group)

long_df %>%
  distinct(PATNO, group_label) %>%
  count(group_label)

# install.packages("ggeffects")   # if not already installed
library(ggeffects)   # for model-based predictions
library(ggplot2)
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
    paste0("`", protein, "` ~ 1 + group_label + time_since_onset + group_label:time_since_onset + (1 + time_since_onset | PATNO)")
  )
  formula_reduced_non_random <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time_since_onset + group_label:time_since_onset + (1 | PATNO)")
  )
  formula_full_slope <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time_since_onset + group_label:time_since_onset + envolved_age_c + (1 + time_since_onset | PATNO)")
  )
  formula_full_non_random <- as.formula(
    paste0("`", protein, "` ~ 1 + group_label + time_since_onset + group_label:time_since_onset + envolved_age_c + (1 | PATNO)")
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
    reason_A <- if (inherits(fit_reduced_slope_ml, "error")) "error_reduced_slope_ML"
    else paste0("reduced_slope_non_random_ML")
  } else {
    fit_A_final_ml <- fit_reduced_slope_ml
    reason_A <- "reduced_slope_ML"
  }
  if (!inherits(fit_A_final_ml, "error") && isSingular(fit_A_final_ml, tol = 1e-4)) {
    reason_A <- paste0(reason_A, "_reduced_non_random_ML")
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
    reason_full <- if (inherits(fit_full_slope_ml, "error")) "error_full_slope_ML"
    else paste0("full_slope_non_random_ML")
  } else {
    fit_full_final_ml <- fit_full_slope_ml
    reason_full <- "full_slope_ML"
  }
  if (!inherits(fit_full_final_ml, "error") && isSingular(fit_full_final_ml, tol = 1e-4)) {
    reason_full <- paste0(reason_full, "_full_non_random_ML")
  }
  
  # 1c. If both fail, skip
  if (inherits(fit_A_final_ml, "error") && inherits(fit_full_final_ml, "error")) {
    warning("Protein '", protein, "' could not be fitted under ML.")
    next
  }
  
  # 1d. Record AICs (for comparison under ML)
  AIC_A_ml    <- if (!inherits(fit_A_final_ml, "error"))    AIC(fit_A_final_ml)    else Inf
  AIC_full_ml <- if (!inherits(fit_full_final_ml, "error")) AIC(fit_full_final_ml) else Inf
  
  if (AIC_full_ml < AIC_A_ml) {
    chosen_model_ml       <- fit_full_final_ml
    chosen_fixed_spec_ml  <- "with_envolved_age_c"
    chosen_random_ml      <- reason_full
  } else {
    chosen_model_ml       <- fit_A_final_ml
    chosen_fixed_spec_ml  <- "without_envolved_age_c"
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
  #best_formula <- formula(chosen_model_ml)
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
    est_time_since_onset_REML            = if("time_since_onset"             %in% rownames(coef_mat)) coef_mat["time_since_onset",             "Estimate"] else NA,
    pval_time_since_onset_REML           = if("time_since_onset"             %in% rownames(coef_mat)) coef_mat["time_since_onset",             "Pr(>|t|)"]  else NA,
    est_group_PD_LRRK2_REML  = if("group_labelPD_LRRK2" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2", "Estimate"] else NA,
    pval_group_PD_LRRK2_REML = if("group_labelPD_LRRK2" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2", "Pr(>|t|)"]  else NA,
    est_group_labelPD_LRRK2_time_since_onset_REML     = if("group_labelPD_LRRK2:time_since_onset" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2:time_since_onset", "Estimate"] else NA,
    pval_group_labelPD_LRRK2_time_since_onset_REML    = if("group_labelPD_LRRK2:time_since_onset" %in% rownames(coef_mat)) coef_mat["group_labelPD_LRRK2:time_since_onset", "Pr(>|t|)"]  else NA,
    est_envolved_age_c_REML    = if("envolved_age_c"     %in% rownames(coef_mat)) coef_mat["envolved_age_c",     "Estimate"] else NA,
    pval_envolved_age_c_REML   = if("envolved_age_c"     %in% rownames(coef_mat)) coef_mat["envolved_age_c",     "Pr(>|t|)"]  else NA,
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
  mutate(
    Protein = sub("^([^_]*_){2}", "", Protein)  # keep text after second underscore
  ) %>%
  # Precompute –log10(pvalue)
  mutate(
    negLog10P = -log10(pvalue)
  )


pdf("The volcano plot of BL_FC and p value with p value.pdf",6,6)
EnhancedVolcano(volcano_df,
                #selectLab = " ",
                lab           = volcano_df$Protein,         
                x             = 'log2FC',    
                y             = 'pvalue',   
                xlab          = bquote(~Log[2]~ "FC"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.1,                    
                FCcutoff      = 0.1,
                xlim = c(-0.5,0.5),
                ylim = c(0,2),
                gridlines.major = F, gridlines.minor = F,
                min.segment.length = 0.0001,
                pointSize = 1.0,
                labSize = 5.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 1,
                legendLabels=c('NS','Log2FC','pvalue','Log2FC & pvalue'),
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 30,
                subtitle = ""
)

dev.off()


volcano_df1 <- df_summaries %>%
  filter(!is.na(est_group_labelPD_LRRK2_time_since_onset_REML) & !is.na(pval_group_labelPD_LRRK2_time_since_onset_REML)) %>%
  select(
    Protein  = protein,
    log2FC   = est_group_labelPD_LRRK2_time_since_onset_REML,
    pvalue   = pval_group_labelPD_LRRK2_time_since_onset_REML
  ) %>%
  mutate(
    Protein = sub("^([^_]*_){2}", "", Protein)  # keep text after second underscore
  ) %>%
  # Precompute –log10(pvalue)
  mutate(
    negLog10P = -log10(pvalue)
  )


pdf("The volcano plot of disease slope_change and p value.pdf",6,6)
EnhancedVolcano(volcano_df1,
                #selectLab = " ",
                lab           = volcano_df1$Protein,              
                x             = 'log2FC',       
                y             = 'pvalue',    
                xlab          = bquote(~"Disease_Slope_Change"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.1,                     
                FCcutoff      = 0.01,
                xlim = c(-0.1,0.1),
                ylim = c(0,2),
                gridlines.major = F, gridlines.minor = F,
                min.segment.length = 0.0001,
                pointSize = 1.0,
                labSize = 5.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 1,
                legendLabels=c('NS','Log2FC','pvalue','Log2FC & pvalue'),
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 15,
                subtitle = ""
)

dev.off()


volcano_df2 <- df_summaries %>%
  filter(!is.na(est_time_since_onset_REML) & !is.na(pval_time_since_onset_REML)) %>%
  select(
    Protein  = protein,
    log2FC   = est_time_since_onset_REML,
    pvalue   = pval_time_since_onset_REML
  ) %>%
  # (Optional) If your protein column has a prefix like "P177_CSF_", remove it:
  mutate(
    Protein = sub("^([^_]*_){2}", "", Protein)  # keep text after second underscore
  ) %>%
  # Precompute –log10(pvalue)
  mutate(
    negLog10P = -log10(pvalue)
  )# handy to precompute


pdf("The volcano plot of HC slope_change and p value.pdf",6,6)
EnhancedVolcano(volcano_df2,
                #selectLab = " ",
                lab           = volcano_df2$Protein,              
                x             = 'log2FC',       
                y             = 'pvalue',    
                xlab          = bquote(~"HC Slope_Change"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.1,                     
                FCcutoff      = 0.01,
                xlim = c(-0.1,0.1),
                ylim = c(0,2),
                gridlines.major = F, gridlines.minor = F,
                min.segment.length = 0.0001,
                pointSize = 1.0,
                labSize = 5.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 1,
                legendLabels=c('NS','Log2FC','pvalue','Log2FC & pvalue'),
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 15,
                subtitle = ""
)

dev.off()





library(ggplot2)
library(lattice)
library(ggeffects)
library(effects)
library(broom)   
library(dplyr)
library(ggplot2)
library(rlang)
library(broom.mixed)

long_df$group_label <- factor(long_df$group_label)


protein_name <- "LMAN2"  #replace the protein you want to plot
proteins  <- paste0("P177_CSF_", protein_name)


# 2. Fit by REML
for (response_var in proteins)  {
  fml <- as.formula(
    paste0(
      " updrs3_score ~ 1 + " , response_var, "*group_label*time_since_onset + envolved_age_c + (1 +", response_var, "|PATNO)"  #+time_since_onset   envolved_age_c +
    )
  )
  m_mixed <- mixed(
    formula = fml,
    data    = long_df, 
    control = lmerControl(optimizer = "bobyqa"),
    method = 'LRT'
  )
  
  m_mixed
  capture.output(
    m_mixed,
    file = paste0(
      response_var,
      "_fix_effect_significance_test_by_mix.txt"
    )
  )
  
  m_reml <- lmer(
    formula = fml,
    data    = long_df, 
    control = lmerControl(optimizer = "bobyqa")
  )
  
  summary(m_reml)
  
  capture.output(
    summary(m_reml),
    file = paste0(
      response_var,
      " _lmer_full_summary.txt"
    )
  )
  
  #  3. Build a prediction grid
  pred_df <- expand.grid(
    time_since_onset  = seq(min(long_df$time_since_onset), max(long_df$time_since_onset), length.out = 100),
    group_label = levels(long_df$group_label)
  )
  # Make sure it’s a factor
  pred_df$group_label <- factor(pred_df$group_label, levels = levels(long_df$group_label))
  
  # 2. Add the centered age and sex columns
  pred_df$envolved_age_c <- 0
  #pred_df$Sex        <- factor("Female", levels = levels(long_df$Sex))
  #str(pred_df)
  
  # 4. Predict the marginal fit
  pred_df$fit <- predict(m_reml, newdata = pred_df, re.form = NA)
  
  # 5. Plot
  # 1. Tidy the model to get estimates + p-values
  tidy_fix <- broom::tidy(m_reml, effects = "fixed")
  
  # extract HC slope (time_since_onset) and its p
  slope_hc   <- tidy_fix  %>% filter(term == "time_since_onset") %>% pull(estimate)
  p_hc       <- tidy_fix  %>% filter(term == "time_since_onset") %>% pull(p.value)
  
  # extract group-specific interaction
  b_int      <- tidy_fix  %>% filter(term == "group_labelPD_LRRK2:time_since_onset") %>% pull(estimate)
  p_int      <- tidy_fix  %>% filter(term == "group_labelPD_LRRK2:time_since_onset") %>% pull(p.value)
  
  # PD_LRRK2 slope = β_time_since_onset + β_group_label:time_since_onset
  slope_pd   <- slope_hc + b_int
  p_pd       <- p_int  # p-value for that interaction
  
  # 2. Build an annotation data frame
  #    We'll place HC text at x= max(time_since_onset)*0.7, y = max(response)*0.9
  #    and PD at x= same, y = max(response)*0.8
  yrange     <- range(long_df[[response_var]], na.rm = TRUE)
  xpos       <- max(long_df$time_since_onset) * 0.7
  annot_df   <- tibble(
    group_label = c("HC", "PD_LRRK2"),
    x           = xpos,
    y           = c(yrange[2]*0.9, yrange[2]*0.8),
    label       = c(
      sprintf("slope=%.3f\np=%.3g",  slope_hc, p_hc),
      sprintf("slope=%.3f\np=%.3g",  slope_pd, p_pd)
    )
  )
  
  # 3. Rebuild your plot, adding geom_text()
  pdf(paste0(response_var, "_ObservedData_Fit.pdf"), 6, 4)
  tryCatch({
    p <- ggplot(long_df,
                aes(time_since_onset, !!sym(response_var), colour = group_label)) +
      geom_point(alpha = 0.3, size = 1) +
      
      geom_text(
        data = annot_df,
        aes(x = x, y = y, label = label, color = group_label),
        hjust = 0, vjust = -2, size = 3.5,
        show.legend = FALSE
      ) +
      
      geom_line(data = pred_df,
                aes(time_since_onset, fit, colour = group_label, group = group_label),
                size = 1) +
      scale_color_manual(values = c(HC = "black", PD_LRRK2 = "red")) +
      labs(
        x     = "Years",
        y     = response_var,
        title = "Observed Data and Estimated Slopes",
        color = "group"
      ) +
      theme_minimal() +
      theme(
        axis.text      = element_text(size = 12),  
        axis.title     = element_text(size = 12),   
        legend.title   = element_text(size = 12),  
        legend.text    = element_text(size = 12)   
      )
    print(p)                  
  }, finally = dev.off())
}  


#plot a set of proteins you are interested in
library(lme4)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(lattice)

long_df$group_label <- factor(long_df$group_label)
protein_names <- c("CTSD","CTSB","CTSL","CTSF","CTSZ","PSAP","TPP1","HEXA","LAMP2","NPC2","HSPA8", "HSPG2",
                   "ORM2","ATRN","SERPING1","FGB","CSF1R", "APOE", "SERPINB6")

for (protein_name in protein_names) {
  
  response_var <- paste0("P177_CSF_", protein_name,"_gc_z")
  
  ## 1 ── fit model ─────────────────────────────────────────────
  fml <- as.formula(
    paste0(response_var,
           paste0("updrs3_score ~ 1 + ", response_var, " * group_label * time_since_onset + envolved_age_c + (1|PATNO)"))
  )
  
  m_reml <- lmer(fml, data = long_df,
                 control = lmerControl(optimizer = "bobyqa"))
  
  ## 2 ── build prediction grid ────────────────────────────────
  pred_df <- expand.grid(
    time_since_onset        = seq(min(long_df$time_since_onset), max(long_df$time_since_onset), length.out = 100),
    group_label = levels(long_df$group_label)
  ) |>
    mutate(
      envolved_age_c = 0,
      Sex          = factor("Female", levels = levels(long_df$Sex)),
      fit          = predict(m_reml, newdata = cur_data(), re.form = NA)
    )
  
  ## 3 ── observed vs fitted plot ──────────────────────────────
  pdf(paste0(response_var, "_ObservedData_Fit.pdf"), 6, 4)
  tryCatch({
    p <- ggplot(long_df,
                aes(time_since_onset, !!sym(response_var), colour = group_label)) +
      geom_point(alpha = .3, size = 1) +
      # 2. group mean at each time_since_onset (line)
      #      stat_summary(
      #        aes(group = group_label, color = group_label),
      #        fun    = mean,
      #        geom   = "line",
      #        size   = 0.5
      #      ) +
      geom_line(data = pred_df,
                aes(time_since_onset, fit, colour = group_label, group = group_label),
                size = 1) +
      scale_color_manual(values = c(HC = "black", PD_LRRK2 = "red")) +
      labs(
        x     = "time_since_onset",
        y     = response_var,
        title = "Observed Data with Group Means and Fitted Trajectories",
        color = "Group"
      ) +
      theme_minimal() +
      theme(
        axis.text      = element_text(size = 16),  
        axis.title     = element_text(size = 18),  
        legend.title   = element_text(size = 18),  
        legend.text    = element_text(size = 16)    
      )
    print(p)                    
  }, finally = dev.off())
  
  ## 4 ── residuals vs fitted ──────────────────────────────────
  ## residuals vs fitted
  pdf(paste0(response_var, "_Residuals_vs_Fitted.pdf"), 6, 6)
  tryCatch({
    print( plot(m_reml) )            # <- explicit print
  }, finally = dev.off())
  
  ## 5 ── QQ plot of random effects ────────────────────────────
  pdf(paste0(response_var, "_QQ_RandomEffects.pdf"), 6, 6)
  tryCatch({
    print( qqmath(ranef(m_reml)) )   # <- explicit print
  }, finally = dev.off())
  
  ## 6 ── ggpredict ribbon (ggsave auto-closes device) ─────────
  p_ci <- ggpredict(m_reml, c("time_since_onset", "group_label")) |> plot()
  ggsave(paste0(response_var, "_group_time_since_onset_CIs.pdf"),
         plot = p_ci, width = 6, height = 6)
  
  ## 7 ── spaghetti + group means ──────────────────────────────
  df_plot <- long_df |>
    select(PATNO, time_since_onset, group_label, value = all_of(response_var))
  
  pdf(paste0(response_var, "_Spaghetti.pdf"), 6, 6)
  tryCatch({
    p2 <- ggplot(df_plot,
                 aes(time_since_onset, value, group = PATNO, colour = group_label)) +
      geom_line(alpha = .25) +
      stat_summary(aes(group = group_label, colour = group_label),
                   fun = mean, geom = "line", size = 1.2) +
      scale_color_manual(values = c(HC = "black", PD_LRRK2 = "red")) +
      theme_minimal()
    print(p2)
  }, finally = dev.off())
}

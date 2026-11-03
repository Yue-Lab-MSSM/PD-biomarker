# =============================================================================
# 04_CSF_LMM_sex_interaction.R
# Longitudinal analysis of CSF proteins in LRRK2 PD and healthy controls by
# linear mixed-effects models with a group x time x sex interaction
# (PPMI Project 177, mass spectrometry)
#
# Manuscript: "Integrated Biofluid Proteomics Identified Dynamic Functional
#              Biomarkers of LRRK2-Linked Parkinson's Disease Progression"
#
# Steps:
#   1. Select HC and LRRK2 PD (PD_LRRK2); time = years since baseline
#      (onset_EVENT_ID); log2(x + 1) transform of protein abundances;
#      age at enrollment centred; sex coded as a factor (Female = reference)
#   2. PCA (after missMDA EM imputation) for quality control
#   3. For each protein, a linear mixed-effects model (lme4/lmerTest, REML):
#        protein ~ group * time * sex + age + (1 + time | participant)
#      P values by Satterthwaite t tests; Benjamini-Hochberg FDR across
#      proteins for each fixed-effect term
#   4. Volcano plots of the baseline group difference (LRRK2 PD vs HC) and
#      the difference in slope (group x time), using nominal P values
#   5. For selected proteins: observed data with model-fitted trajectories,
#      model diagnostics (residuals vs fitted, Q-Q plot of random effects),
#      predicted trajectories with 95% CIs, and spaghetti plots
#
# Input : data/Working_file_ppmi_project_177_csf.csv
#           (one row per sample: phenotypes and P177_CSF_* protein columns)
# Output: per-protein model details; fixed effects, random effects and
#         FDR-corrected results for all proteins; volcano plots;
#         trajectory and diagnostic plots
# =============================================================================

library(lme4)
library(tidyverse)
library(afex)
library(dplyr)
# Example data preparation (replace with your actual data loading)
df <- read.csv("Working_file_ppmi_project_177_csf without ST_with years after onset.csv")

df$time <- df$onset_EVENT_ID
# Check for unmapped EVENT_ID values
if (any(is.na(df$time))) {
  warning("Some EVENT_ID values could not be mapped to time. Check your data.")
}

# View the first few rows to confirm
head(df[, c("PATNO", "onset_EVENT_ID", "time")])

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
long_df <- df %>% filter(group_label %in% c("HC", "PD_LRRK2")) #for LRRK2_NMC, change the name here

# Identify protein columns (adjust the pattern based on your column names)
protein_cols <- grep("^P177_CSF_", names(long_df), value = TRUE)

# Apply log2 transformation to all protein columns
long_df <- long_df %>%
  mutate(across(all_of(protein_cols), ~ log2(. + 1), .names = "{.col}"))

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
summary(long_df$ENROLL_AGE)

# 3. Mean‐center ENROLL_AGE
long_df$ENROLL_AGE_c <- long_df$ENROLL_AGE - mean(long_df$ENROLL_AGE, na.rm = TRUE)

# 4. Re‐code Sex as factor
#    assuming SexM1.F0 is coded 1 = male, 0 = female
long_df$Sex <- factor(long_df$SexM1.F0,
                      levels = c(0, 1),
                      labels = c("Female", "Male"))


####Generate the values for all the proteins together.
library(lme4)
library(lmerTest)
library(dplyr)

# Initialize a list to store results for each protein
results <- list()
model_details <- list()

# === 1. Fit models & collect results ===
for (protein in protein_cols) {
  # Build the model formula
  formula_str <- paste0(
    protein,
    " ~ 1 + group_label + time + Sex",
    " + group_label:time + time:Sex + group_label:Sex",
    " + group_label:time:Sex + ENROLL_AGE_c",
    " + (1 + time | PATNO)"
  )
  model <- lmer(
    as.formula(formula_str),
    data    = long_df,
    control = lmerControl(optimizer = "bobyqa")
  )
  sm <- summary(model)
  
  # Helper to safely extract estimate/p-value
  getEst <- function(term) tryCatch(sm$coefficients[term, "Estimate"], error = function(e) NA)
  getP   <- function(term) tryCatch(sm$coefficients[term, "Pr(>|t|)"], error = function(e) NA)
  
  # Store key fixed effects + REML
  results[[protein]] <- tibble(
    protein                  = protein,
    coef_group_HC            = getEst("group_labelHC"),
    p_group_HC               = getP("group_labelHC"),
    coef_group_PD_LRRK2      = getEst("group_labelPD_LRRK2"),
    p_group_PD_LRRK2         = getP("group_labelPD_LRRK2"),
    coef_time                = getEst("time"),
    p_time                   = getP("time"),
    coef_SexMale             = getEst("SexMale"),
    p_SexMale                = getP("SexMale"),
    coef_gl_time             = getEst("group_labelPD_LRRK2:time"),
    p_gl_time                = getP("group_labelPD_LRRK2:time"),
    coef_time_Sex            = getEst("time:SexMale"),
    p_time_Sex               = getP("time:SexMale"),
    coef_gl_Sex              = getEst("group_labelPD_LRRK2:SexMale"),
    p_gl_Sex                 = getP("group_labelPD_LRRK2:SexMale"),
    coef_three_way           = getEst("group_labelPD_LRRK2:time:SexMale"),
    p_three_way              = getP("group_labelPD_LRRK2:time:SexMale"),
    coef_ENROLL_AGE_c        = getEst("ENROLL_AGE_c"),
    p_ENROLL_AGE_c           = getP("ENROLL_AGE_c"),
    REML                     = sm$devcomp$cmp["REML"]
  )
  
  # Scaled residuals
  res_qs <- quantile(residuals(model), probs = c(0, .25, .5, .75, 1))
  # Random effects table
  rand_eff <- as.data.frame(VarCorr(model)) %>%
    transmute(
      grp      = grp,
      term     = if_else(is.na(var1), "(Intercept)", var1),
      variance = vcov,
      stddev   = sdcor
    )
  
  # Save detailed pieces
  model_details[[protein]] <- list(
    protein          = protein,
    reml_criterion   = sm$devcomp$cmp["REML"],
    n_obs            = sm$devcomp$dims["n"],
    n_groups         = sm$ngrps[["PATNO"]],
    scaled_residuals = res_qs,
    random_effects   = rand_eff,
    fixed_effects    = as.data.frame(sm$coefficients) %>%
      rownames_to_column("term") %>%
      rename(
        Estimate  = "Estimate",
        Std_Error = "Std. Error",
        DF        = "df",
        T_Value   = "t value",
        P_Value   = "Pr(>|t|)"
      )
  )
}

# === 2. Assemble & FDR-correct ===
results_df <- bind_rows(results) %>%
  mutate(
    fdr_group_HC       = p.adjust(p_group_HC,             method = "fdr"),
    fdr_group_PD_LRRK2 = p.adjust(p_group_PD_LRRK2,      method = "fdr"),
    fdr_time           = p.adjust(p_time,                method = "fdr"),
    fdr_Sex            = p.adjust(p_SexMale,             method = "fdr"),
    fdr_gl_time        = p.adjust(p_gl_time,             method = "fdr"),
    fdr_time_Sex       = p.adjust(p_time_Sex,            method = "fdr"),
    fdr_gl_Sex         = p.adjust(p_gl_Sex,              method = "fdr"),
    fdr_three_way      = p.adjust(p_three_way,           method = "fdr"),
    fdr_ENROLL_AGE     = p.adjust(p_ENROLL_AGE_c,        method = "fdr")
  )

# Identify any proteins with at least one FDR < .05
significant_proteins <- filter(
  results_df,
  if_any(starts_with("fdr_"), ~ . < 0.05)
)

# === 3. Write CSV outputs ===
# 3.1 Summary results
write.csv(results_df,             "fdr_results_all_proteins.csv",  row.names = FALSE)
write.csv(significant_proteins,   "significant_proteins.csv",      row.names = FALSE)

# 3.2 Per-protein detailed CSVs
for (protein in names(model_details)) {
  det <- model_details[[protein]]
  
  # A. Model summary
  summary_df <- tibble(
    Protein          = det$protein,
    REML_Criterion   = det$reml_criterion,
    Num_Observations = det$n_obs,
    Num_Groups       = det$n_groups
  )
  # B. Residuals
  resid_df <- as.data.frame(t(det$scaled_residuals)) %>%
    setNames(paste0("Scaled_Residual_", names(det$scaled_residuals)))
  # C. Random effects
  rand_df <- det$random_effects
  # D. Fixed effects
  fix_df  <- det$fixed_effects
  
  # Combine into one multi‐section CSV
  sink(paste0("model_details_", protein, ".csv"))
  cat("=== Model_Summary ===\n")
  write.csv(summary_df, row.names = FALSE)
  cat("\n=== Scaled_Residuals ===\n")
  write.csv(resid_df,    row.names = FALSE)
  cat("\n=== Random_Effects ===\n")
  write.csv(rand_df,     row.names = FALSE)
  cat("\n=== Fixed_Effects ===\n")
  write.csv(fix_df,      row.names = FALSE)
  sink()
}

# 3.3 Combined details across all proteins
all_summaries   <- bind_rows(lapply(model_details, function(x) {
  tibble(
    Protein          = x$protein,
    REML_Criterion   = x$reml_criterion,
    Num_Observations = x$n_obs,
    Num_Groups       = x$n_groups
  )
}))
all_residuals   <- bind_rows(lapply(model_details, function(x) {
  as.data.frame(t(x$scaled_residuals)) %>%
    setNames(paste0("Scaled_Residual_", names(x$scaled_residuals))) %>%
    mutate(Protein = x$protein)
}))
all_random      <- bind_rows(lapply(model_details, function(x) {
  mutate(x$random_effects, Protein = x$protein)
}))
all_fixed       <- bind_rows(lapply(model_details, function(x) {
  mutate(x$fixed_effects,  Protein = x$protein)
}))

write.csv(all_summaries, "all_proteins_model_summary.csv",   row.names = FALSE)
write.csv(all_residuals, "all_proteins_scaled_residuals.csv", row.names = FALSE)
write.csv(all_random,    "all_proteins_random_effects.csv",   row.names = FALSE)
write.csv(all_fixed,     "all_proteins_fixed_effects.csv",    row.names = FALSE)


##Volcano Plot of Effect Sizes
library(EnhancedVolcano)

library(dplyr)

# fixed_effects_all has columns: Protein, term, estimate, p_value, etc.
volcano_df <- all_fixed %>%
  filter(term == "group_labelPD_LRRK2") %>%      # keep only LRRK2 effect
  select(Protein, log2FC = estimate, pvalue = p_value) %>%
  mutate(Protein = sub("^([^_]*_){2}", "", Protein),  #only get the text after the second "_"
         negLog10P = -log10(pvalue))              # handy to precompute


pdf("The volcano plot of BL_FC and p value with p value.pdf",6,6)
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
                # blank‐out the default labels:
                #title    = "",      # no main title
                subtitle = "",      # no subtitle
                #caption  = " "       # no caption (if you want)
)

dev.off()



# fixed_effects_all has columns: Protein, term, estimate, p_value, etc.
volcano_df1 <- all_fixed %>%
  filter(term == "group_labelPD_LRRK2:time") %>%      # keep only LRRK2 effect
  select(Protein, log2FC = estimate, pvalue = p_value) %>%
  mutate(Protein = sub("^([^_]*_){2}", "", Protein),  #only get the text after the second "_"
         negLog10P = -log10(pvalue))              # handy to precompute


pdf("The volcano plot of slop_change and p value.pdf",6,6)
EnhancedVolcano(volcano_df1,
                #selectLab = " ",
                lab           = volcano_df1$Protein,               # name of the column with labels
                x             = 'log2FC',       # the log₂FC column
                y             = 'pvalue',    # the p-value column
                xlab          = bquote(~"Slop_Change"),
                ylab          = bquote(~-Log[10]~italic(P)),
                pCutoff       = 0.1,                     # raw p-value cutoff
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
                max.overlaps = 30,
                # blank‐out the default labels:
                #title    = "",      # no main title
                subtitle = "",      # no subtitle
                #caption  = " "       # no caption (if you want)
)

dev.off()

library(ggplot2)
library(lattice)
library(ggeffects)
library(effects)
# 1. Ensure Group is a factor

long_df$group_label <- factor(long_df$group_label)


protein_name <- "CTSD" #change to the proteins you are interested in
response_var  <- paste0("P177_CSF_", protein_name)

# 2. Fit by REML

# build a one-off formula string and convert to a formula object

fml <- as.formula(
  paste0(
    response_var,
    " ~ 1 + group_label + time + Sex + group_label:time +  
    time:Sex + group_label:Sex + group_label:time:Sex + ENROLL_AGE_c + (1+time|PATNO)"
  )
)

m_reml <- lmer(
  formula = fml,
  data    = long_df, 
  control = lmerControl(optimizer = "bobyqa")
)

summary(m_reml)

#  3. Build a prediction grid
pred_df <- expand.grid(
  time  = seq(min(long_df$time), max(long_df$time), length.out = 100),
  group_label = levels(long_df$group_label)
)
# Make sure it’s a factor
pred_df$group_label <- factor(pred_df$group_label, levels = levels(long_df$group_label))

# 2. Add the centered age and sex columns
pred_df$ENROLL_AGE_c <- 0
pred_df$Sex        <- factor("Female", levels = levels(long_df$Sex))
#str(pred_df)

# 4. Predict the marginal fit
pred_df$fit <- predict(m_reml, newdata = pred_df, re.form = NA)


# 5. Plot

pdf(
  paste0(response_var ,"_Observed Data with Fitted Trajectories.pdf"),6,6
)

ggplot(long_df, aes(x = time, y = !!sym(response_var), color = group_label)) +
  geom_point(alpha = 0.3, size = 1) +
  
  # 2. group mean at each time (line)
  #  stat_summary(
  #    aes(group = group_label, color = group_label),
  #    fun    = mean,
  #    geom   = "line",
  #    size   = 0.5
  #  )
  
  # 3. fitted marginal trajectories
  geom_line(
    data = pred_df,
    aes(x = time, y = fit, group = group_label),
    size = 1
  ) +
  
  geom_line(data = pred_df, aes(y = fit), size = 1) +
  labs(
    x = "Time",
    y = response_var,
    title = "Observed Data with Fitted Trajectories"
  ) +
  
  scale_color_manual(
    values = c(
      HC        = "black",    # or whatever palette you like
      PD_LRRK2  = "red"
    )
  ) +
  theme_minimal()
dev.off()

###Diagnostic Plots for Model Checking
#Residuals vs. Fitted and Q–Q Plots.  Purpose: Check homoscedasticity, linearity, and normality of residuals.

pdf(
  paste0(response_var ,"_residuals vs. fitted.pdf"),6,6
)
plot(m_reml)            # base R in lmer: residuals vs. fitted
dev.off()


pdf(
  paste0(response_var ,"_lattice Q–Q for random effects.pdf"),6,6
)
qqmath(ranef(m_reml))   # lattice Q–Q for random effects
dev.off()


pdf(
  paste0(response_var ,"_group_label_time_CIs.pdf"),6,6
)

ggpredict(m_reml, c("time", "group_label")) %>% plot()
dev.off()

pdf(
  paste0(response_var ,"_seperate_group_label_time_CIs.pdf"),6,6
)

plot(effect("group_label:time", m_reml))
dev.off()


pdf(
  paste0("The protein plot of", response_var ,".pdf"),6,6
)

#Trajectory (“Spaghetti”) Plots of Individual Proteins
prot <- response_var
df <- long_df %>%
  select(PATNO, time, group_label, value = all_of(prot))

ggplot(df, aes(x = time, y = value, group = PATNO, color = group_label)) +
  # thin, transparent lines for each subject
  geom_line(alpha = 0.3) +
  #ylim(16,21) +
  
  # thick lines for the group means, inheriting the same colour mapping
  stat_summary(
    aes(group = group_label, color = group_label), 
    fun    = mean, 
    geom   = "line", 
    size   = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      HC        = "black",    # or whatever palette you like
      PD_LRRK2  = "red"
    )
  ) +
  
  labs(
    title    = paste("Time Course of", prot),
    x        = "Time",
    y        = "Log₂(abundance)",
    color    = "Group"
  ) +
  theme_minimal()

dev.off()



#plot a set of protiens
library(lme4)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(lattice)

long_df$group_label <- factor(long_df$group_label)
protein_names <- c("CTSD","ORM2","LYZ","SERPING1",
                   "CTSB","CTSL","CTSF","CTSZ")

for (protein_name in protein_names) {
  
  response_var <- paste0("P177_CSF_", protein_name)
  
  ## 1 ── fit model ─────────────────────────────────────────────
  fml <- as.formula(
    paste0(response_var,
           " ~ 1 + group_label + time + Sex + group_label:time +  
    time:Sex + group_label:Sex + group_label:time:Sex + ENROLL_AGE_c + (1+time|PATNO)")
  )
  
  m_reml <- lmer(fml, data = long_df,
                 control = lmerControl(optimizer = "bobyqa"))
  
  ## 2 ── build prediction grid ────────────────────────────────
  pred_df <- expand.grid(
    time        = seq(min(long_df$time), max(long_df$time), length.out = 100),
    group_label = levels(long_df$group_label)
  ) |>
    mutate(
      ENROLL_AGE_c = 0,
      Sex          = factor("Female", levels = levels(long_df$Sex)),
      fit          = predict(m_reml, newdata = cur_data(), re.form = NA)
    )
  
  ## 3 ── observed vs fitted plot ──────────────────────────────
  pdf(paste0(response_var, "_ObservedData_Fit.pdf"), 6, 6)
  tryCatch({
    p <- ggplot(long_df,
                aes(time, !!sym(response_var), colour = group_label)) +
      geom_point(alpha = .3, size = 1) +
      # 2. group mean at each time (line)
      stat_summary(
        aes(group = group_label, color = group_label),
        fun    = mean,
        geom   = "line",
        size   = 0.5
      ) +
      geom_line(data = pred_df,
                aes(time, fit, colour = group_label, group = group_label),
                size = 1) +
      scale_color_manual(values = c(HC = "black", PD_LRRK2 = "red")) +
      labs(
        x     = "Time",
        y     = response_var,
        title = "Observed Data with Group Means and Fitted Trajectories",
        color = "Group"
      ) +
      theme_minimal()
    print(p)                      # <- required inside pdf()
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
  p_ci <- ggpredict(m_reml, c("time", "group_label")) |> plot()
  ggsave(paste0(response_var, "_group_time_CIs.pdf"),
         plot = p_ci, width = 6, height = 6)
  
  ## 7 ── spaghetti + group means ──────────────────────────────
  df_plot <- long_df |>
    select(PATNO, time, group_label, value = all_of(response_var))
  
  pdf(paste0(response_var, "_Spaghetti.pdf"), 6, 6)
  tryCatch({
    p2 <- ggplot(df_plot,
                 aes(time, value, group = PATNO, colour = group_label)) +
      geom_line(alpha = .25) +
      stat_summary(aes(group = group_label, colour = group_label),
                   fun = mean, geom = "line", size = 1.2) +
      scale_color_manual(values = c(HC = "black", PD_LRRK2 = "red")) +
      theme_minimal()
    print(p2)
  }, finally = dev.off())
}

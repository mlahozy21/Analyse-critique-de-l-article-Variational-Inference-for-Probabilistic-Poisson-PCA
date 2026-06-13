# ==============================================================================
# MASTER SCRIPT: REPRODUCTION OF SECTION 7.2 (OAK PATHOBIOME)
# Dataset: Oaks (Jakuschkin et al., 2016)
# Reference: Chiquet et al. (2018), Figures 4 & 5
# ==============================================================================

# 1. Load Libraries
library(PLNmodels)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

message("--- Step 1: Data Loading & Preparation ---")

# Load data
data(oaks)

# --- CRITICAL FIX: Prepare Covariates Manually ---
# The covariates in the 'oaks' object are loose vectors. 
# We must group them into a data frame.
my_covariates <- data.frame(
  tree = oaks$tree,
  orientation = oaks$orientation,
  distTOground = oaks$distTOground
)

# --- CRITICAL FIX: Assign Row Names ---
# We must assign the row names from the abundance matrix to the covariates
# so that PLNmodels knows they belong to the same samples.
rownames(my_covariates) <- rownames(oaks$Abundance)

# Create the PLN data object (including Offsets for normalization)
oaks_data <- prepare_data(
  counts = oaks$Abundance, 
  covariates = my_covariates, 
  offset = oaks$Offset
)

message("--- Step 2: Fitting Models (Ranks 1 to 30) ---")

# Fit Model M0: Null model (Offset only, no covariates)
# Reference: "We first fitted a null Poisson-lognormal model M0" 
# Capture failures explicitly: a try({}) with no error handler would let the
# script continue with a stale/undefined M0_collection and silently mask the
# failure. tryCatch re-raises with a clear message instead.
M0_collection <- tryCatch(
  PLNPCA(Abundance ~ 1, data = oaks_data, ranks = 1:30),
  error = function(e) stop("PLNPCA fit for M0 failed: ", conditionMessage(e))
)

# Fit Model M1: Full model (Offset + Tree + Orientation)
# Reference: "Model M1 involving two covariates: tree... and orientation" 
M1_collection <- tryCatch(
  PLNPCA(Abundance ~ tree + orientation, data = oaks_data, ranks = 1:30),
  error = function(e) stop("PLNPCA fit for M1 failed: ", conditionMessage(e))
)

message("Models fitted successfully. Generating figures...")


# ==============================================================================
# FIGURE 4: MODEL SELECTION & DIAGNOSTICS (Complete)
# Matches 
# ==============================================================================

# Helper function to safely extract R2 or calculate it if missing
get_criteria_robust <- function(model_collection, label) {
  df <- model_collection$criteria
  df$model_type <- label
  df$k <- 1:nrow(df)
  
  # Robust extraction of R2 (handles cases where column is named differently or missing)
  r2_values <- sapply(model_collection$models, function(x) {
    val <- x$R_squared
    if(is.null(val)) val <- x$R2
    if(is.null(val)) val <- NA 
    return(val)
  })
  df$R_squared <- as.numeric(r2_values)

  # If R2 is genuinely unavailable from the fitted models we must NOT fabricate
  # a substitute. The previous "fallback" min-max normalised the log-likelihood
  # and relabelled it "R2" — that quantity is not an R2 (it is a rescaled
  # log-lik that is 0 at the worst rank and 1 at the best by construction), and
  # plotting it on an R2 axis is misleading. Error out honestly instead so the
  # problem is fixed at the source (e.g. enable rsquared in the post-treatment
  # config) rather than silently masked.
  if (all(is.na(df$R_squared))) {
    stop("R2 is missing for every model in '", label, "'. Refit with the ",
         "post-treatment R2 computation enabled (config_post: rsquared = TRUE) ",
         "instead of substituting a fabricated criterion.")
  }

  return(df)
}

# Extract data
df_M0 <- get_criteria_robust(M0_collection, "M0")
df_M1 <- get_criteria_robust(M1_collection, "M1")
all_criteria <- rbind(df_M0, df_M1)

# Calculate Entropy Approximation: BIC - ICL 
all_criteria$Entropy <- all_criteria$BIC - all_criteria$ICL

# --- FIGURE 4a: Model Selection Criteria (Left: M0, Right: M1) ---

# Prepare data for M0 plot
df_M0_long <- df_M0 %>%
  select(k, loglik, BIC, ICL) %>%
  pivot_longer(cols = c("loglik", "BIC", "ICL"), names_to = "Criterion", values_to = "Value")

# Plot Left (M0)
p4a_left <- ggplot(df_M0_long, aes(x = k, y = Value, color = Criterion)) +
  geom_line(linewidth = 0.8) + geom_point() + theme_bw() +
  labs(title = "Offset (M0)", 
       x = "Number of axes", y = "Criterion value") +
  theme(legend.position = "none") # No legend to avoid duplication

# Prepare data for M1 plot
df_M1_long <- df_M1 %>%
  select(k, loglik, BIC, ICL) %>%
  pivot_longer(cols = c("loglik", "BIC", "ICL"), names_to = "Criterion", values_to = "Value")

# Plot Right (M1)
p4a_right <- ggplot(df_M1_long, aes(x = k, y = Value, color = Criterion)) +
  geom_line(linewidth = 0.8) + geom_point() + theme_bw() +
  labs(title = "Offset and covariates (M1)", 
       x = "Number of axes", y = "") +
  theme(legend.position = c(0.7, 0.3))

# --- FIGURE 4b: Goodness of Fit and Entropy ---

# Left: R2 Criterion
p4b_left <- ggplot(all_criteria, aes(x = k, y = R_squared, color = model_type)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) + geom_point(na.rm = TRUE) + theme_bw() + 
  ylim(0, 1) +
  labs(title = "R2 Criterion", 
       x = "Number of axes", y = "R2") +
  theme(legend.position = "none")

# Right: Entropy
p4b_right <- ggplot(all_criteria, aes(x = k, y = Entropy, color = model_type)) +
  geom_line(linewidth = 0.8) + geom_point() + theme_bw() +
  labs(title = "Entropy (BIC - ICL)", 
       x = "Number of axes", y = "Entropy") +
  theme(legend.position = "bottom")

# Assemble Figure 4 (Top row: 4a, Bottom row: 4b)
grid.arrange(p4a_left, p4a_right, p4b_left, p4b_right, ncol = 2, 
             top = "Figure 4: Model Selection & Diagnostics")


# ==============================================================================
# FIGURE 5: INDIVIDUAL FACTOR MAPS (PCA)
# Comparison of M0 (uncorrected) vs M1 (corrected for tree effect)
# Matches 
# ==============================================================================

# Select best models based on ICL 
best_M0 <- getBestModel(M0_collection, crit = "ICL")
best_M1 <- getBestModel(M1_collection, crit = "ICL")

message(paste("Selected Rank for M0:", best_M0$rank))
message(paste("Selected Rank for M1:", best_M1$rank))

# --- Top Row: Effect of 'Tree' ---

# M0 (Null) - Colored by Tree
p5_top_left <- plot(best_M0, map = "individual", ind_cols = oaks$tree) +
  ggtitle("Offset (M0)", subtitle = "Color: Tree Status") + 
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Axis 1", y = "Axis 2")

# M1 (Covariates) - Colored by Tree
p5_top_right <- plot(best_M1, map = "individual", ind_cols = oaks$tree) +
  ggtitle("Offset + Covariates (M1)", subtitle = "Color: Tree Status") + 
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Axis 1", y = "Axis 2")


# --- Bottom Row: Effect of 'Distance to Ground' ---

# M0 (Null) - Colored by Distance
p5_bot_left <- plot(best_M0, map = "individual", ind_cols = oaks$distTOground) +
  ggtitle("Offset (M0)", subtitle = "Color: Dist. to Ground") + 
  scale_color_viridis_c() + 
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Axis 1", y = "Axis 2")

# M1 (Covariates) - Colored by Distance
p5_bot_right <- plot(best_M1, map = "individual", ind_cols = oaks$distTOground) +
  ggtitle("Offset + Covariates (M1)", subtitle = "Color: Dist. to Ground") + 
  scale_color_viridis_c() + 
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Axis 1", y = "Axis 2")

# Assemble Figure 5
grid.arrange(p5_top_left, p5_top_right, p5_bot_left, p5_bot_right, ncol = 2, 
             top = "Figure 5: PCA Projections - Uncorrected (M0) vs Corrected (M1)")

message("Done! Figures 4 and 5 generated.")

# ==============================================================================
# REPRODUCTION OF FIGURE 6: THE "NUCLEAR OPTION" (Unique Variable Names)
# ==============================================================================

library(PLNmodels)
library(ggplot2)
library(dplyr)

message("--- Step 1: Clean Slate Data Prep ---")

# Load fresh data
data(oaks)

# Unique name for covariates dataframe
covs_final <- data.frame(
  tree = oaks$tree,
  orientation = oaks$orientation,
  distTOground = oaks$distTOground
)
rownames(covs_final) <- rownames(oaks$Abundance)

# Unique name for PLN data object
data_final <- prepare_data(
  counts = oaks$Abundance, 
  covariates = covs_final, 
  offset = oaks$Offset
)

message("--- Step 2: Fitting Rank 25 (Unique Object) ---")

# We use a totally NEW variable name 'collection_final'
# This prevents R from looking at 'M1_optimal' or 'M1_collection' from before
collection_final <- tryCatch(
  PLNPCA(Abundance ~ tree + orientation, data = data_final, ranks = 25),
  error = function(e) stop("PLNPCA fit (rank 25, Figure 6) failed: ", conditionMessage(e))
)

# Access the first (and only) model in this new collection
# We call it 'model_25' to be specific
model_25 <- collection_final$models[[1]]

message("--- Step 3: Extracting variational parameters ---")

# Conditional standard error of the latent position (see 02_conditional_variance.R
# for the full derivation). The variational posterior of the latent score vector
# of sample i is Z_i ~ N(M_i, diag(S_i^2)); the latent position in species space
# is Z C^T with C = model_par$C the (p x q) PCA loadings. Hence the conditional
# variance of the (i, j) latent position is ( S^2 %*% t(C^2) )_{ij}, an (n x p)
# matrix aligned one-to-one with the (n x p) abundance matrix Y.
#
# FIX: the previous code used B = model_par$B (the covariate-regression
# coefficients, d x p) instead of the PCA loadings C, and then transposed /
# truncated matrices to "common rank" to force compatible shapes. B plays no
# role in the latent-score covariance; the correct quantity uses S and C only,
# and the shapes are correct by construction (no reshaping hacks needed).
C_final <- as.matrix(model_25$model_par$C)   # (p x q) loadings, NOT model_par$B
S_final <- as.matrix(model_25$var_par$S)     # (n x q) variational std-devs
Y_final <- as.matrix(data_final$Abundance)   # (n x p) observed abundances

cat("dim C (p x q):", paste(dim(C_final), collapse = " x "), "\n")
cat("dim S (n x q):", paste(dim(S_final), collapse = " x "), "\n")

# Shapes must be correct by construction; assert rather than silently reshape.
stopifnot(ncol(C_final) == 25L, ncol(S_final) == 25L)
stopifnot(nrow(C_final) == ncol(Y_final))    # p = number of species
stopifnot(nrow(S_final) == nrow(Y_final))    # n = number of samples


# --- Step 4: Calculation & Plot ---

# Var[(Z C^T)_ij | Y] = ( S^2 %*% t(C^2) )_ij   ->   SE = sqrt(.)
S_sq <- S_final^2
C_sq <- C_final^2

Var_Z <- S_sq %*% t(C_sq)   # (n x p), aligns with Y
SD_Z  <- sqrt(Var_Z)

# DataFrame for ggplot
df_plot <- data.frame(
  Abundance = as.vector(Y_final),
  SD_Z      = as.vector(SD_Z)
)

p6_final <- ggplot(df_plot, aes(x = Abundance, y = SD_Z)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_x_log10(breaks = c(1, 10, 100, 500, 1000, 2000),
                labels = c("1", "10", "100", "500", "1k", "2k")) +
  theme_bw() +
  labs(
    title = "Figure 6: Conditional Standard Error of Z",
    subtitle = "Model M1 (Rank 25) - Final Reproduction",
    x = "Abundance (log scale)",
    y = "Approximate conditional standard error"
  ) +
  annotation_logticks(sides = "b")

print(p6_final)
message("Figure 6 is finally generated.")
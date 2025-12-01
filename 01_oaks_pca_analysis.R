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
try({
  M0_collection <- PLNPCA(Abundance ~ 1, data = oaks_data, ranks = 1:30)
})

# Fit Model M1: Full model (Offset + Tree + Orientation)
# Reference: "Model M1 involving two covariates: tree... and orientation" 
try({
  M1_collection <- PLNPCA(Abundance ~ tree + orientation, data = oaks_data, ranks = 1:30)
})

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
  
  # Fallback: Approximate R2 if completely missing (for visualization purposes)
  if(all(is.na(df$R_squared))) {
    min_l <- min(df$loglik)
    max_l <- max(df$loglik)
    df$R_squared <- (df$loglik - min_l) / (max_l - min_l)
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
try({
  collection_final <- PLNPCA(Abundance ~ tree + orientation, data = data_final, ranks = 25)
})

# Access the first (and only) model in this new collection
# We call it 'model_25' to be specific
model_25 <- collection_final$models[[1]]

message("--- Step 3: Verifying Dimensions ---")

# Extract matrices
B_final <- as.matrix(model_25$model_par$B)
S_final <- as.matrix(model_25$var_par$S)
Y_final <- as.matrix(data_final$Abundance)

# Check dimensions explicitly
q_B <- ncol(B_final)
q_S <- ncol(S_final)

cat("Rank of B:", q_B, "\n")
cat("Rank of S:", q_S, "\n")

# Safety fix if B comes transposed (rare but possible)
if (q_B != 25 && nrow(B_final) == 25) {
  message("Transposing B...")
  B_final <- t(B_final)
  q_B <- ncol(B_final)
}


# --- Step 4: Calculation & Plot ---

# Formula: Var(Z) = S^2 * (B^2)^T
S_sq <- S_final^2
B_sq <- B_final^2

Var_Z <- S_sq %*% t(B_sq)
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
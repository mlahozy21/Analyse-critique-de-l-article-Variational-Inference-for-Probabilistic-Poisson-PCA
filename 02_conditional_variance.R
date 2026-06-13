# ==============================================================================
# REPRODUCTION OF FIGURE 6: EXACT RANK (q=25)
# Reference: Section 7.2 states optimal rank for M1 is 25.
#
# Figure 6 (Chiquet et al., 2018) plots, for every (sample i, species j) cell,
# the *conditional standard error* of the latent position against the observed
# abundance Y_ij. In PLNPCA the variational posterior of the latent score vector
# of sample i is Gaussian:  Z_i ~ N(M_i, diag(S_i^2)),  with M, S of size (n x q)
# (var_par$M and var_par$S; S holds the variational *standard deviations*, S2
# their squares). The latent position in the p-dimensional species space is
# Z C^T, where C = model_par$C are the (p x q) PCA loadings. Because the
# variational covariance of Z_i is diagonal, the conditional variance of the
# (i, j) latent position is
#
#       Var[(Z C^T)_ij | Y] = sum_k C_{jk}^2 S_{ik}^2  =  ( S^2 %*% t(C^2) )_{ij}
#
# i.e. an (n x p) matrix whose square root is the per-cell conditional SE that
# aligns one-to-one with the (n x p) abundance matrix Y. (This is exactly the
# 0.5 * S2 %*% t(C2) variance-correction term that appears in the package's
# fitted intensity A = exp(Z + 0.5 * S2 %*% t(C2)); see PLNPCAfit-class.R.)
#
# PREVIOUS BUG: the earlier version computed  S^2 %*% t(B^2)  using B =
# model_par$B, the *covariate-regression* coefficients (d x p), NOT the PCA
# loadings C (p x q). That is dimensionally and conceptually wrong: B has
# nothing to do with the latent-score covariance, and the shapes only "lined
# up" because of transpose/truncate-to-common-rank reshaping hacks that forced
# the matrices into compatible sizes. Those hacks have been removed; the
# correct quantity is a pure function of S and C.
# ==============================================================================

# 1. CLEAN SLATE
rm(list = ls())

library(PLNmodels)
library(ggplot2)
library(dplyr)

message("--- Step 1: Loading & Preparing Data ---")
data(oaks)

# Fix row names mapping
covs <- data.frame(
  tree = oaks$tree,
  orientation = oaks$orientation,
  distTOground = oaks$distTOground
)
rownames(covs) <- rownames(oaks$Abundance)

oaks_data <- prepare_data(
  counts = oaks$Abundance,
  covariates = covs,
  offset = oaks$Offset
)

message("--- Step 2: Forcing Rank 25 ---")

# We fit specifically Rank 25. Capture failures explicitly so we never proceed
# on a stale / non-existent model object.
model_collection <- tryCatch(
  PLNPCA(Abundance ~ tree + orientation, data = oaks_data, ranks = 25),
  error = function(e) stop("PLNPCA fit (rank 25) failed: ", conditionMessage(e))
)

# Access the generated model directly (it is the only one in the list)
my_model <- model_collection$models[[1]]

message("--- Step 3: Extracting variational parameters ---")

# C : (p x q) PCA loadings        (model_par$C)
# S : (n x q) variational std-dev (var_par$S; S^2 = var_par$S2 is the variance)
# Y : (n x p) observed abundances
C <- as.matrix(my_model$model_par$C)   # loadings, NOT model_par$B
S <- as.matrix(my_model$var_par$S)
Y <- as.matrix(oaks_data$Abundance)

# Sanity checks on the shapes (no silent transpose/truncation: the dimensions
# must be correct by construction, otherwise something upstream is wrong).
stopifnot(ncol(C) == 25L, ncol(S) == 25L)          # both have q = 25 columns
stopifnot(nrow(C) == ncol(Y))                      # C is (p x q), p = #species
stopifnot(nrow(S) == nrow(Y))                      # S is (n x q), n = #samples

# 4. Conditional standard error per (sample, species) cell.
#    Var[(Z C^T)_ij | Y] = ( S^2 %*% t(C^2) )_ij   -> SE = sqrt(.)
S_sq <- S^2          # (n x q) variational variances
C_sq <- C^2          # (p x q) squared loadings

Var_Z <- S_sq %*% t(C_sq)   # (n x p), aligns with Y
SD_Z  <- sqrt(Var_Z)

# 5. Plotting (Exact Paper Style)
df_plot <- data.frame(
  Abundance = as.vector(Y),
  SD_Z      = as.vector(SD_Z)
)

# We use coord_cartesian to match the paper's zoom (0 to 0.8 approx)
p6 <- ggplot(df_plot, aes(x = Abundance, y = SD_Z)) +
  geom_point(size = 0.6, alpha = 0.4) +
  scale_x_log10(breaks = c(1, 10, 100, 500, 1000, 2000),
                labels = c("1", "10", "100", "500", "1k", "2k")) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.0)) + # Zoom in to match Figure 6 limits
  labs(
    title = paste("Figure 6: Conditional Standard Error"),
    subtitle = "Matches Chiquet et al. (2018) structure",
    x = "Abundance Y (log scale)",
    y = "Conditional standard error of the latent position"
  ) +
  annotation_logticks(sides = "b")

print(p6)
message("SUCCESS: Figure 6 generated")

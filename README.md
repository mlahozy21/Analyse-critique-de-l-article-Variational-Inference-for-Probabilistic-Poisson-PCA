# Reproduction of "Variational Inference for Probabilistic Poisson PCA" (Chiquet et al., 2018)

This repository contains the **R code** required to reproduce the experimental results presented in the paper:

> **Chiquet, J., Mariadassou, M., & Robin, S. (2018).** Variational inference for probabilistic Poisson PCA. *The Annals of Applied Statistics*, 12(4), 2674-2698.

## Overview

The objective of this project is to validate the **Probabilistic Poisson PCA (PLNPCA)** method using the official R package `PLNmodels`. The scripts focus on **Section 7.2 (Oak Powdery Mildew Pathobiome)**, reproducing the model selection diagnostics, PCA projections, and conditional variance analysis.

The code addresses several implementation challenges, including:
* Manual preparation of covariate data frames to resolve package version incompatibilities.
* Robust extraction of model criteria ($R^2$, BIC, ICL).
* Manual matrix algebra to compute conditional standard errors for Figure 6.

## Repository Structure

| File Name | Description | Corresponds to |
| :--- | :--- | :--- |
| **`01_oaks_pca_analysis.R`** | Main analysis script. Fits Null Model ($M_0$) vs. Covariate Model ($M_1$), performs model selection, and visualizes the latent space. | **Figures 4 & 5** |
| **`02_conditional_variance.R`** | Advanced analysis script. Fits the optimal rank model ($q=25$) and computes the approximate conditional standard error of latent variables $Z_{ij}$. | **Figure 6** |

## 🛠️ Prerequisites

To run these scripts, you need **R** installed along with the following packages:

```r
install.packages(c("PLNmodels", "ggplot2", "dplyr", "tidyr", "gridExtra"))

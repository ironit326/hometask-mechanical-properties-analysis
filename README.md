# Statistical Analysis of AM50 Magnesium Alloy Mechanical Properties

## Overview
This R script performs comprehensive statistical analysis of tensile test data for AM50 magnesium alloy, comparing Reference and Nb-B inoculated materials with different sample geometries (Flat 5mm and Round 6mm).

## Data Description
- **YS** - Yield Strength (MPa)
- **UTS** - Ultimate Tensile Strength (MPa)  
- **A** - Elongation at break (%)
- **Unknown** - Unknown parameter
- **Type** - Material type (Ref / Nb-B)
- **Size** - Sample size (5mm / 6mm)
- **Metrics** - Sample geometry (Flat / Round)

## Analyses Performed
1. **Descriptive statistics** (mean, median, SD, skewness, kurtosis)
2. **Normality testing** (Shapiro-Wilk test)
3. **Histograms** with optimized bin count (Sturges' rule)
4. **Boxplots** for group comparison
5. **Scatter plots** with regression lines
6. **ANOVA** (two-way)
7. **Correlation matrix** (Pearson)
8. **Linear regression** with VIF, R², RMSE

## Required Libraries
```r
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(corrplot)
library(car)
library(caret)

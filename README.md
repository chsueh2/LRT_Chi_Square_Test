# Chi Square Test for Homogeneity: Likelihood Ratio Test vs. Monte Carlo Simulation

Use Likelihood Ratio Test (LRT) to conduct a chi-square test for homogeneity and compare it with a Monte Carlo simulation.

[Project report](https://rpubs.com/clh2021/1113681)

Key features:

- Likelihood Ratio Test
- Multinomial Distribution
- Chi-square Test for Homogeneity
- Pearson Chi-square statistic
- $\alpha$ Control
- Asymptotic Test Power
- Monte Carlo Simulation 

R packages used:

- `conflicted`: tools to deal with R package conflicts
- `here`: enables easy file referencing and builds file paths in a OS-independent way
- `scales`: formats and labels scales nicely for better visualization
- `skimr`: provides summary statistics about variables in data frames, tibbles, data tables and vectors
- `glue`: embeds and evaluates R expressions into strings to be printed as messages
- `gt`: pretty-looking tables
- `GGally`: extends 'ggplot2' by combining geometric objects with transformed data
- `moments`: calculate moments, kurtosis and skewness
- `cumstats`: efficiently compute cumulative standard deviation
- `tidyverse`: includes collections of useful packages like `dplyr` (data manipulation), `tidyr` (tidying data),  `ggplots` (creating graphs), etc.
- `zealot`: provides a `%<-%` operator to perform multiple, unpacking, and destructing assignment 

## Project Report

[Project report](https://rpubs.com/clh2021/1113681) ([Github Markdown](./project3.md))([PDF](./project3.pdf))([R Markdown](./project3.Rmd)).

The analysis results with all theoretical backgrounds and math derivations are included.

Chien-Lan Hsueh (chienlan.hsueh at gmail.com)

## Overview and Project Goal

Consider three different hospitals and each hospital has patients that end up with infections. We are interested whether or not the multinomials are homogenous across the hospitals.

Conduct chi-square tests for homogeneity using R. Also manually calculate the LRT statistic, the Pearson Chi-square statistic, the critical value, and find approximate p-values for the hypotheses using both test statistics.

- Determine how well the asymptotic rejection region performs at controlling $\alpha$
- Determine the power of the asymptotic test when comparing certain alternative situations

## Part 1 - Data Example

Create matrix for the given hospital data and conduct chi-square tests for homogeneity. Implement the  test in R and verify the result with R build-in test function `chisq.test()`.

## Part 2 - Derivation

Derive the likelihood ratio test by considering the generic case of comparing $J$ independent multinomial distributions, each with $I$ categories. 

## Part 3 - Simulation

The Pearson statistic can be derived as a Taylor series approximation to the LRT. Investigate the $\alpha$ control of the Pearson chi-square test and its power. 

Consider the following three different setups:

1. Two multinomial case only, where each multinomial has three categories
1. All combinations of four sample sizes for each multinomial (16 total cases)
1. Three different probabilities that may generate a particular multinomial

## Part 4 - Summarization 

Analyze and summarize the results using faceted line plots for easy comparison.


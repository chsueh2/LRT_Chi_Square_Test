# Chi Square Test for Homogeneity: Likelihood Ratio Test vs. Monte Carlo Simulation


# Packages ----------------------------------------------------------------

if (!require("pacman")) utils::install.packages("pacman", dependencies = TRUE)

pacman::p_load(
  conflicted, here, 
  scales, skimr, glue, gt, GGally, 
  moments, cumstats,
  tidyverse, zeallot
)

# resolve conflicts by setting preference
conflict_prefer("select", "dplyr", quiet = T)
conflict_prefer("filter", "dplyr", quiet = T)

# infix operator: string concatenation (ex: "act" %&% "5")
'%&%' <- function(x, y) paste0(x, y)

# infix operator: not %in%
'%notin%' <- Negate('%in%')


# Helper Functions --------------------------------------------------------

# calculate expected counts from a contingency table
expected_counts <- function(obs){
  # convert a vector to a nx1 matrix
  obs <- as.matrix(obs)
  
  # column sums (as 1xm matrix) and row sums (as nx1 matrix)
  csum <- colSums(obs) %>% as.matrix() %>% t()
  rsum <- rowSums(obs) %>% as.matrix()
  
  # expected counts
  exp <- rsum %*% csum / sum(csum)
  # add 0.5 to any expected counts that end up being 0
  exp[exp == 0] <- 0.5
  
  return(exp)
}

# LRT or Pearson test statistic
chisq_stat <- function(obs, exp = expected_counts(obs), type = "Pearson"){
  # return Pearson's test stat by default
  switch (type,
          LRT = list(type = "LRT", test_stat = 2*sum(obs*log(obs/exp))),
          list(type = "Pearson", test_stat = sum((obs-exp)^2/exp))
  )
}

# chi-square test
chisq_test_GOF <- function(obs, alpha = 0.05, type = "Pearson"){
  # degree of freedom
  ## df = k-1 if one-dimension contingency table
  ## df = (R-1)(C-1) if two-dimension contingency table
  df <- prod(dim(as.matrix(obs))-1)
  df <- ifelse(df == 0, length(obs)-1, df)
  
  # get test stat
  c(type, test_stat) %<-% chisq_stat(obs, type = type)
  
  # p-value and critical value
  p_value <- pchisq(test_stat, df, lower.tail = F)
  critical <- qchisq(alpha, df, lower.tail = F)
  
  # print out the test summary
  print(glue(
    "{type} Chi-squared test\n",
    "X-squared = {round(test_stat, 3)}, df = {df},",
    "p-value = {label_scientific(4)(p_value)}\n",
    "alpha = {alpha}, critical value = {round(critical, 3)}\n",
    "Reject H0? {p_value < alpha}"
  ))
  
  # return results
  list(
    type = type,
    test_stat = test_stat,
    df = df,
    p_value = p_value,
    alpha = alpha,
    crtical = critical)
}

# simulate data and get empirical rejection rate
rejection_rate <- function(N, n1, p1, n2, p2, cutoff, type = "Pearson", seed = 2022){
  # make the simulation reproducible 
  set.seed(seed)
  
  # simulate N contingency tables from 2 multinomials
  data <- replicate(
    N, 
    cbind(rmultinom(1, n1, p1), rmultinom(1, n2, p2)) %>% t(), 
    simplify = F)
  
  # return the empirical rejection rate
  ## `chisq_stat(x, type = type)$test_stat` gives Pearson's or LRT test statistics
  mean(sapply(data, function(x){ chisq_stat(x, type = type)$test_stat }) > cutoff, na.rm = T)
}

# Data --------------------------------------------------------------------

# create a matrix for the example data
example_obs <- matrix(
  data = c(41, 27, 51, 36, 3, 40, 169, 106, 109),
  nrow = 3, 
  byrow = TRUE,
  dimnames = list(
    LETTERS[1:3],
    c("Surgical Site Infections", "Pneumonia Infections", "Bloodstream Infections")
  )
)

example_obs

# perform Pearson's chi-square test
test_Pearson <- chisq_test_GOF(example_obs)

# double check our results with R's `chisq.test()`
chisq.test(example_obs)

# perform LRT chi-square test
test_LRT <- chisq_test_GOF(example_obs, type = "LRT")


# Simulation --------------------------------------------------------------

# seed for random number generator
seed <- 2022

# significance level and cutoff
alpha <- 0.05
cutoff <- qchisq(alpha, df = 2, lower.tail = F)

# simulation size: # of contingency tables in each study case
N <- 50*1000

# sample sizes: total counts in each contingency table
sample_sizes <- c(20, 30, 50, 100)

# probabilities of multinomial distribution
probs <- list(c(1, 1, 1), c(1, 3, 6), c(1, 1, 8))

# types of chi-square test statistics
types <- c("Pearson", "LRT")

# MC simulation for empirical alpha
df_alpha <- expand_grid(
  n1 = sample_sizes,
  p1 = probs,
  n2 = sample_sizes,
  type = types
) %>% 
  rowwise() %>% 
  mutate(
    alpha = rejection_rate(N, n1, p1, n2, p2 = p1, cutoff, type, seed),
    # create legend labels
    match = match(list(p1), probs),
    prob = c("Equal", "Mixed1", "Mixed2")[match]
  )

# print the results
df_alpha %>% 
  select(-match) %>% mutate(across(everything(), toString))

# plot comparison
df_alpha %>% 
  ggplot(aes(n1, alpha, col = prob)) +
  geom_line(linewidth = 2) +
  facet_grid(rows = vars(type), cols = vars(n2)) +
  labs(
    title = "Empirical alpha",
    subtitle = glue("Significance level: {alpha}")
  )

# MC simulation for empirical alpha
df_power <- expand_grid(
  p1 = probs,
  n1 = sample_sizes,
  n2 = sample_sizes,
  p2 = probs,
  type = types
) %>% 
  rowwise() %>% 
  mutate(
    match1 = match(list(p1), probs),
    match2 = match(list(p2), probs)
  ) %>% 
  # remove unwanted combinations
  filter(match1 < match2) %>% 
  mutate(
    power = rejection_rate(N, n1, p1, n2, p2, cutoff, type, seed),
    # create legend labels
    prob = glue("{match1} vs {match2}") 
  )

# print the results
df_power %>% 
  select(-match1, -match2) %>% 
  mutate(across(everything(), toString))

# plot comparison
df_power %>% 
  ggplot(aes(n1, power, col = prob)) +
  geom_line(linewidth = 2) +
  facet_grid(rows = vars(type), cols = vars(n2)) +
  labs(
    title = "Empirical Power",
    subtitle = glue("Significance level: {alpha}")
  )


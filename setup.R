# source("miscfunctions.R")
source("helperfunctions.R")
library(lme4)
library(cowplot)
library(rstan)
library(rstanarm)
library(egg)
evidscheme <- c("Baseline"="black","Inculpatory"="#ba4040ff","Exculpatory"="#406bbaff", "Ambiguous"="#765884ff")
theme_set(cowplot::theme_minimal_grid(font_size=11))

shrinky_beta <- stan_model("models/shrinkaget_beta.stan")
shrinky_beta_div <- stan_model("models/shrinkaget_beta_divnorm.stan")
shrinky_beta_divreg <- stan_model("models/shrinkaget_beta_divreg.stan")
shrinky_beta_divreg_byval <- stan_model("models/shrinkaget_beta_divreg_wbyval.stan")

shrinky <- stan_model("models/shrinkaget.stan")
shrinky_div <- stan_model("models/shrinkaget_divnorm.stan")
shrinky_divreg <- stan_model("models/shrinkaget_divreg.stan")
shrinky_divreg_byval <- stan_model("models/shrinkaget_divreg_wbyval.stan")
shrinky_divreg_contrast <- stan_model("models/shrinkaget_divreg_contrast_coded.stan")

bdat <- readindat("exculpatory_burdenofproof", pthresh=0.05)
bdat_cond <- readindat("exculpatory_conditional", pthresh=0.05)[cond_evidence %in% c("balanced", "credible")]
bdat_rate <- readindat("exculpatory_rateless", pthresh=0.05)[cond_evidence %in% c("balanced", "credible")]

evnames <- model.matrix(~ physical + document + witness + character, data=bdat)[,-1] |> colnames()

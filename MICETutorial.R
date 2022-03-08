## ----setup, include=FALSE--------------------------------------------------------
library(tidyverse); library(broom);
library(knitr); library(ggplot2);
library(geomtextpath);
library(mice);
library(miceadds);
library(survival);
knitr::opts_chunk$set(echo = T, warning = F, message = F);
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size);
})

def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize", paste0("\n \\", options$size,"\n\n", x, "\n\n \\normalsize"), x)
})

#knitr::opts_chunk$set(width = 10);
#knitr::opts_chunk$set(tidy.opts=list(width.cutoff=40));
recache = F;
options(digits = 4);
figure_scaler = 1/2; #1/2 for ioslides; ~1/3 for word, pdf
text_scaler = 3/3;#1 for ioslides; 2/3 for word, pdf
fig.x = 19 * figure_scaler;
fig.y = 10.5 * figure_scaler;


## ---- out.width = "95%",echo = F-------------------------------------------------
knitr::include_graphics("MiceAlgorithm.png")


## ---- out.width = "95%",echo = F-------------------------------------------------
knitr::include_graphics("NormalImputations.png")


## ---- out.width = "80%",echo = F-------------------------------------------------
knitr::include_graphics("donors.png")


## --------------------------------------------------------------------------------
args(parlmice)


## ---- include = T, echo = T, cache = T-------------------------------------------
library(survival)
cgd_first_obs <- 
  cgd %>% 
  # Keep only first event
  filter(enum == 1) %>%
  mutate(bmi = weight / (height / 100) ^2,
         nih = str_detect(center, "NIH")) %>%
  select(tstop, status, treat, sex, age, weight, height, bmi, inherit, propylac, nih)


## ---- include = T, echo = T, cache = T-------------------------------------------
head(cgd_first_obs)


## ---- include = T, echo = T, cache = T-------------------------------------------

ggplot(cgd_first_obs) + 
  geom_point(aes(x = height, y = weight, color = (bmi < 80)))

ggplot(cgd_first_obs) + 
  geom_point(aes(x = age, y = weight, color = (bmi < 80)))

ggplot(cgd_first_obs) + 
  geom_point(aes(x = age, y = height, color = (bmi < 80)))

coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac, data = cgd_first_obs %>% filter(bmi < 80))

coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac, data = cgd_first_obs)



## ---- include = T, echo = T, cache = T-------------------------------------------
# set to missing the strange value of height
cgd_first_obs <- 
  cgd_first_obs %>% 
  mutate(height = 
           case_when(bmi < 80 ~ height, 
                     TRUE ~ NA_real_),
         bmi = weight / (height / 100) ^2)


## ----full_model, include = T, echo = T, cache = T--------------------------------
# no missingness (except that one strange height)
full_model <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac, data = cgd_first_obs)

full_coefs <- 
  tidy(full_model, conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "full")

full_coefs %>% 
  mutate(across(where(is.numeric), ~formatC(.x, format = "f", digits = 3))) %>%
  kable()


## ---- include = T, echo = T, cache = T-------------------------------------------

permute_true_false <- 
  function(n, prob) {
    if(is.logical(n)) {
      # n can be a vector of trues and falses
      result = rep(F, length(n))
      result[n] = sample(c(rep(T, round(prob*sum(n))),rep(F, sum(n) - round(prob*sum(n))) ))
    } else if(length(n) == 1 && n%%1 == 0) {
      # otherwise n should be a single integer
      result = sample(c(rep(T, round(prob*n)),rep(F, n - round(prob*n)) ))
    } else {stop("'n' must be a vector of logicals or a single integer")}
    result
  }

set.seed(1)
cgd_first_obs_mcar <- 
  cgd_first_obs %>% 
  mutate(
    # If TRUE, then keep it; otherwise it will be missing
    height = case_when(
      permute_true_false(n(), 0.8) ~ height),
    inherit = case_when(
      permute_true_false(n(), 0.75) ~ inherit),
    treat = case_when(
      permute_true_false(n(), 0.75) ~ treat),
    bmi = weight / (height / 100) ^2)

md.pattern(cgd_first_obs_mcar)

# do a dry run of mice to prefill and then edit mice arguments
init = mice(cgd_first_obs_mcar, m = 1, max = 0, print = F);
method = init$meth;
method;
# Formulas for derived variables;
method["bmi"] = 
  "~I(weight / (height / 100) ^2)";
method;

predictorMatrix = init$pred;

predictorMatrix

# Don't use bmi (a derived variable) to predict height (its ingredient):
predictorMatrix[ "height", "bmi"] = 0;
# Don't use ingredients of derived variables to predict anything else:
predictorMatrix[,c("height", "weight")] = 0;
predictorMatrix["height", "weight"] = 1;
predictorMatrix["weight", "height"] = 1;
# Ensure derived variables depend only on ingredients
predictorMatrix["bmi",] = 0;
predictorMatrix["bmi", c("weight", "height")] = 1;

predictorMatrix

visitSequence = init$visitSequence;

visitSequence

visitSequence = c(c("height", "bmi"),setdiff(visitSequence, c("height", "bmi")));

visitSequence

n_mi_chain_per_core = 2;#1;#
n_cores = 5; # should be < then value of parallel::detectCores()
n_mi_chain = n_mi_chain_per_core * n_cores; 
# Total number of imputed datasets will be n_mi_chain_per_core * n_cores
n_iter =  30;#10;#

cgd_first_obs_mcar_imputed =  
  parlmice(cgd_first_obs_mcar,
           predictorMatrix = predictorMatrix,
           method = method,
           visitSequence = visitSequence,
           maxit = n_iter,
           printFlag = T, 
           n.imp.core = n_mi_chain_per_core,
           cluster.seed = 1,
           n.core = n_cores, 
           cl.type = "FORK")

mcar_mi_coefs <- 
  with(cgd_first_obs_mcar_imputed, coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac)) %>%
  pool() %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "mcar_mi")

mcar_si_coefs <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac,
        data = complete(cgd_first_obs_mcar_imputed, 1)) %>%
  tidy(conf.int = T) %>%
  select(term,  estimate, conf.low, conf.high) %>%
  mutate(approach = "mcar_si")

mcar_cc_coefs <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac,
        data = cgd_first_obs_mcar) %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "mcar_cc")



## ---- include = T, echo = T, cache = T-------------------------------------------
# Stack all imputed datasets
complete(cgd_first_obs_mcar_imputed, "long", include = T) %>% 
  # add the actual dataset for comparison
  bind_rows(cgd_first_obs %>% select(colnames(.)) %>% 
              mutate(.imp = -1)) %>%
  # group by type of data
  mutate(.imp_category = 
           case_when(.imp == -1 ~ "1full",
                     .imp == 0  ~ "2obs", 
                     .imp > 0 ~ "3imp")) %>%
  group_by(.imp_category) %>% 
  summarize(median_bmi = median(bmi, na.rm = T), 
            upperq_bmi = quantile(bmi, p = 0.75, na.rm = T),
            meanrIFNg = mean(treat == "rIFN-g", na.rm = T), 
            mean_autosomal= mean(inherit == "autosomal", na.rm = T))


## ---- echo = F, cache = F, size = "scriptsize", fig.width = fig.x, fig.height=fig.y----
bind_rows(full_coefs,mcar_cc_coefs, mcar_si_coefs, mcar_mi_coefs) %>%
  mutate(approach = factor(approach) %>% fct_inorder()) %>%
  mutate(ci_width = conf.high - conf.low) %>%
  ggplot(aes(x = estimate, y = approach)) + 
  geom_point() + 
  geom_label(aes(label = formatC(ci_width, format = "f", digits = 2), 
                 x = estimate + (0.8) * ci_width)) + 
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  facet_wrap(~term, scales = "free")



## ---- include = T, echo = T, cache = T-------------------------------------------

set.seed(1)
cgd_first_obs_mar <- 
  cgd_first_obs %>% 
  mutate(
    # This is stil MCAR as before
    height = case_when(
      permute_true_false(n(), 0.8) ~ height),
    # This is MAR as long as we condition on sex in the imputation
    inherit = case_when(
      sex == "female" & permute_true_false(sex == "female", 0.4) ~ inherit,
      sex == "male" & permute_true_false(sex == "male", 0.83) ~ inherit),
    # This is MAR as long as we condition on nih in the imputation
    treat = case_when(
      nih & permute_true_false(nih, 0.16) ~ treat,
      !nih & permute_true_false(!nih, 0.90)  ~ treat),
    bmi = weight / (height / 100) ^2)

cgd_first_obs_mar_imputed =  
  parlmice(cgd_first_obs_mar,
           predictorMatrix = predictorMatrix,
           method = method,
           visitSequence = visitSequence,
           maxit = n_iter,
           printFlag = T, 
           n.imp.core = n_mi_chain_per_core,
           cluster.seed = 1,
           n.core = n_cores, 
           cl.type = "FORK")


## ---- include = T, echo = F, cache = T-------------------------------------------

mar_mi_coefs <- 
  with(cgd_first_obs_mar_imputed, coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac)) %>%
  pool() %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "mar_mi")

mar_si_coefs <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac,
        data = complete(cgd_first_obs_mar_imputed, 1)) %>%
  tidy(conf.int = T) %>%
  select(term,  estimate, conf.low, conf.high) %>%
  mutate(approach = "mar_si")



## ---- include = T, echo = T, cache = T-------------------------------------------

mar_cc_coefs <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac,
        data = cgd_first_obs_mar) %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "mar_cc")



## ---- include = T, echo = T, cache = T-------------------------------------------

complete(cgd_first_obs_mar_imputed, "long", include = T) %>% 
  bind_rows(cgd_first_obs %>% select(colnames(.)) %>% 
              mutate(.imp = -1)) %>%
  mutate(.imp_category = 
           case_when(.imp == -1 ~ "1full",
                     .imp == 0  ~ "2obs", 
                     .imp > 0 ~ "3imp")) %>%
  group_by(.imp_category) %>% 
  summarize(median_bmi = median(bmi, na.rm = T), 
            upperq_bmi = quantile(bmi, p = 0.75, na.rm = T),
            meanrIFNg = mean(treat == "rIFN-g", na.rm = T), 
            mean_autosomal= mean(inherit == "autosomal", na.rm = T))



## ---- echo = F, cache = F, size = "scriptsize", fig.width = fig.x, fig.height=fig.y----
bind_rows(full_coefs,mar_cc_coefs, mar_si_coefs, mar_mi_coefs) %>%
  mutate(approach = factor(approach) %>% fct_inorder()) %>%
  mutate(ci_width = conf.high - conf.low) %>%
  ggplot(aes(x = estimate, y = approach)) + 
  geom_point() + 
  geom_label(aes(label = formatC(ci_width, format = "f", digits = 2), 
                 x = estimate + (0.8) * ci_width)) + 
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  facet_wrap(~term, scales = "free")



## ---- include = T, echo = T, cache = T-------------------------------------------
set.seed(1)
cgd_first_obs_nmar <- 
  cgd_first_obs %>% 
  mutate(
    # This is stil MCAR as before
    height = case_when(
      permute_true_false(n(), 0.8) ~ height),
    # This is NMAR: the probability of missingness depends on the missing value
    inherit = case_when(
      inherit == "autosomal" & permute_true_false(inherit == "autosomal", 0.95) ~ inherit,
      inherit == "X-linked" & permute_true_false(inherit == "X-linked", 0.65) ~ inherit),
    # This is the same as above but will be NMAR when we don't include
    # nih in the imputation step
    treat = case_when(
      nih & permute_true_false(nih, 0.16) ~ treat,
      !nih & permute_true_false(!nih, 0.90)  ~ treat),
    bmi = weight / (height / 100) ^2)

# If we don't adjust for nih in our imputation models, then 
# the mar scenario becomes an nmar scenario:
predictorMatrix_nmar = predictorMatrix
predictorMatrix_nmar[,"nih"] = 0

cgd_first_obs_nmar_imputed =  
  parlmice(cgd_first_obs_nmar,
           predictorMatrix = predictorMatrix_nmar,
           method = method,
           visitSequence = visitSequence,
           maxit = n_iter,
           printFlag = T, 
           n.imp.core = n_mi_chain_per_core,
           cluster.seed = 1,
           n.core = n_cores, 
           cl.type = "FORK")

# Compare to scenario where we *do* adjust for nih
cgd_first_obs_mar_imputed =  
  parlmice(cgd_first_obs_nmar,
           predictorMatrix = predictorMatrix,
           method = method,
           visitSequence = visitSequence,
           maxit = n_iter,
           printFlag = T, 
           n.imp.core = n_mi_chain_per_core,
           cluster.seed = 1,
           n.core = n_cores, 
           cl.type = "FORK")


## ---- include = T, echo = F, cache = T-------------------------------------------

mar_mi_coefs <- 
  with(cgd_first_obs_mar_imputed, coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac)) %>%
  pool() %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "mar_mi")

nmar_mi_coefs <- 
  with(cgd_first_obs_nmar_imputed, coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac)) %>%
  pool() %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "nmar_mi")

nmar_si_coefs <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac,
        data = complete(cgd_first_obs_nmar_imputed, 1)) %>%
  tidy(conf.int = T) %>%
  select(term,  estimate, conf.low, conf.high) %>%
  mutate(approach = "nmar_si")


nmar_cc_coefs <- 
  coxph(Surv(tstop, status) ~ treat + sex + age + bmi + inherit + propylac,
        data = cgd_first_obs_nmar) %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(approach = "nmar_cc")



## ---- include = T, echo = T, cache = T-------------------------------------------
complete(cgd_first_obs_nmar_imputed, "long", include = T) %>% 
  bind_rows(cgd_first_obs %>% select(colnames(.)) %>% 
              mutate(.imp = -1)) %>%
  mutate(.imp_category = 
           case_when(.imp == -1 ~ "1full",
                     .imp == 0  ~ "2obs", 
                     .imp > 0 ~ "3imp")) %>%
  group_by(.imp_category) %>% 
  summarize(median_bmi = median(bmi, na.rm = T), 
            upperq_bmi = quantile(bmi, p = 0.75, na.rm = T),
            meanrIFNg = mean(treat == "rIFN-g", na.rm = T), 
            mean_autosomal= mean(inherit == "autosomal", na.rm = T))


complete(cgd_first_obs_mar_imputed, "long", include = T) %>% 
  bind_rows(cgd_first_obs %>% select(colnames(.)) %>% 
              mutate(.imp = -1)) %>%
  mutate(.imp_category = 
           case_when(.imp == -1 ~ "1full",
                     .imp == 0  ~ "2obs", 
                     .imp > 0 ~ "3imp")) %>%
  group_by(.imp_category) %>% 
  summarize(median_bmi = median(bmi, na.rm = T), 
            upperq_bmi = quantile(bmi, p = 0.75, na.rm = T),
            meanrIFNg = mean(treat == "rIFN-g", na.rm = T), 
            mean_autosomal= mean(inherit == "autosomal", na.rm = T))



## ---- echo = F, cache = F, size = "scriptsize", fig.width = fig.x, fig.height=fig.y----
bind_rows(full_coefs,nmar_cc_coefs, nmar_si_coefs, nmar_mi_coefs, mar_mi_coefs) %>%
  mutate(approach = factor(approach) %>% fct_inorder()) %>%
  mutate(ci_width = conf.high - conf.low) %>%
  ggplot(aes(x = estimate, y = approach)) + 
  geom_point() + 
  geom_label(aes(label = formatC(ci_width, format = "f", digits = 2), 
                 x = estimate + (0.8) * ci_width)) + 
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  facet_wrap(~term, scales = "free")



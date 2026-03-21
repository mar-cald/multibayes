# worked example

# load libraries
library(afex); library(brms); library(dplyr)

# Select data from experiment 1; see ? afex::stroop ------
data_stroop <- subset(afex::stroop, subset = study == 1)|> 
  # remove useless levels id
  dplyr::mutate(id = droplevels(pno)) |> 
  # select relevant variables
  dplyr::select(id, condition, congruency, acc, rt)|>
  # remove NA 
  na.omit()



## Reaction time analysis--------------------------
data_rt <- data_stroop |>
  group_by(id) |>
  # remove error and post error trials
  filter(acc == 1 & lag(acc, default = 1) == 1) |>
  ungroup()

formula <- rt ~ congruency*condition + 
  (1 + congruency*condition || id)

mod_rt <- brm(formula = formula, 
           family = lognormal(),
           cores = 3, chains = 3,
           iter = 10000,
           save_pars = save_pars(all = FALSE, group = FALSE),
           data = data_rt)

save(mod_rt, file = "paper/script/example/mod_rt.rda")


## Accuracy analysis--------------------------------

data_acc_agg <- data_stroop |>
  group_by(id, condition, congruency) |>
  summarise(hits = sum(acc), n = n(), .groups = "drop")

formula <- hits | trials(n) ~ 1 + congruency*condition + 
  (1 + congruency*condition || id)

mod_acc <- brm(formula = formula, 
           family = binomial(link = "logit"),
           cores = 3, chains = 3,
           iter = 10000, 
           data = data_acc_agg)

save(mod_acc, file = "paper/script/example/mod_acc.rda")



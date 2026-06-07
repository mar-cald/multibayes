# worked example

# load libraries
library(afex); library(brms); library(dplyr)

# Select data from experiment 1; see ? afex::stroop ------
data_stroop <- subset(afex::stroop, subset = study == 1 & 
                        pno %in% sample(levels(pno), size = 80))|> 
  # remove useless levels id
  dplyr::mutate(id = droplevels(pno)) |> 
  # select relevant variables
  dplyr::select(id, condition, congruency, trialnum, acc, rt)|>
  # previous cong
  group_by(id, condition)|>
  mutate(prev_cong = lag(congruency),
         trial = scale(trialnum)[,1])|>
  ungroup()|>
  na.omit()


## Reaction time analysis--------------------------
data_rt <- data_stroop |>
  group_by(id) |>
  # remove error and post error trials
  filter(acc == 1 & lag(acc, default = 1) == 1) |>
  ungroup()

contrasts(data_rt$congruency) <- contr.sum(2)/2
contrasts(data_rt$prev_cong) <- contr.sum(2)/2
contrasts(data_rt$condition) <- contr.sum(2)/2

formula <- rt ~ 1 + congruency * prev_cong + 
  congruency * trial * condition +
  (1 + congruency + prev_cong + 
     congruency + trial + condition || id)

mod_rt <- brm(formula = formula, 
           family = lognormal(),
           cores = 3, chains = 3,
           iter = 8000,
           save_pars = save_pars(all = FALSE, group = FALSE),
           data = data_rt)

save(mod_rt, file = "paper/script/example/mod_rt.rda")



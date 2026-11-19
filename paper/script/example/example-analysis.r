# worked example
# run from the repository root: Rscript paper/script/example/example-analysis.r

# load libraries
library(afex); library(brms); library(dplyr)

# all 253 participants of Study 1; see ? afex::stroop ------
stroop1 <- subset(afex::stroop, subset = study == 1)

data_stroop <- stroop1 |>
  # remove useless levels id
  dplyr::mutate(id = droplevels(pno)) |>
  # select relevant variables
  dplyr::select(id, condition, congruency, trialnum, acc, rt) |>
  # previous trial within each block (participant x condition),
  # computed before any trial is removed
  arrange(id, condition, trialnum) |>
  group_by(id, condition) |>
  mutate(prev_cong = lag(congruency), prev_acc = lag(acc)) |>
  ungroup()


## Reaction time analysis--------------------------
# correct trials preceded by a correct trial; this also removes trials without
# a response (acc = NA), trials after a missing response, and the first trial
# of each block
data_rt <- data_stroop |>
  filter(acc == 1, prev_acc == 1)

contrasts(data_rt$congruency) <- contr.sum(2)/2
contrasts(data_rt$prev_cong) <- contr.sum(2)/2
contrasts(data_rt$condition) <- contr.sum(2)/2
# trial standardized (z-score)
data_rt$trial <- as.numeric(scale(data_rt$trialnum))

formula <- rt ~ 1 + congruency * prev_cong +
  congruency * condition * trial +
  (1 + prev_cong + congruency * condition + trial * condition || id)

mod_rt <- brm(formula = formula,
              family = lognormal(),
              cores = 3, chains = 3,
              iter = 8000, seed = 1,
              save_pars = save_pars(all = FALSE, group = FALSE),
              data = data_rt)

save(mod_rt, file = "paper/script/example/mod_rt.rda")

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

formula <- rt ~ condition*congruency + 
  (1 + condition*congruency || id)

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

formula <- hits | trials(n) ~ 1 + condition*congruency + 
  (1 + condition*congruency || id)

mod_acc <- brm(formula = formula, 
           family = binomial(link = "logit"),
           cores = 3, chains = 3,
           iter = 10000, 
           data = data_acc_agg)

save(mod_acc, file = "paper/script/example/mod_acc.rda")


# Extract and combine draws ------------------------
# Extract draws
draws_rt  <- as.data.frame(mod_rt)
draws_acc <- as.data.frame(mod_acc)

# check colnames
head(colnames(draws_rt), n = 4)
par_names_rt <- colnames(draws_rt)[2:4]

head(colnames(draws_acc), n = 4)
par_names_acc <- colnames(draws_acc)[2:4]

# Select fixed-effect columns of interest
draws_rt  <- draws_rt[,  par_names_rt]
draws_acc <- draws_acc[, par_names_acc]

# Combine draws across models
draws <- cbind(draws_rt, draws_acc)
colnames(draws) <- c("condition_rt", "congruency_rt", 
                     "condition:congruency_rt", "condition_acc", "congruency_acc", 
                     "condition:congruency_acc")

pd.adjust(draws = draws, q = 0.4, R = TRUE)

cor(draws_acc)[lower.tri(cor(draws_acc))]
cor(draws_rt)[lower.tri(cor(draws_rt))]
cor(draws)[lower.tri(cor(draws))]

df = cor(draws)

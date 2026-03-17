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
  

# condition (low demand = -0.5; high demand = 0.5) 
# Stroop congruency (congruent = -0.5; incongruent = 0.5). 
contrasts(data_stroop$condition) <- -(contr.sum(2)/2)
contrasts(data_stroop$congruency) <- -(contr.sum(2)/2)


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
par_names <- colnames(draws_rt)[2:4]

# Select fixed-effect columns of interest
draws_rt  <- draws_rt[,  par_names]
draws_acc <- draws_acc[, par_names]

# Combine draws across models
draws <- cbind(draws_rt, draws_acc)
colnames(draws) <- c("condition_rt", "congruency_rt", 
                     "condition:congruency_rt", "condition_acc", "congruency_acc", 
                     "condition:congruency_acc")

pd.adjust(draws = draws, q = 0.4, R = TRUE)




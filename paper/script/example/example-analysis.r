# worked example

# load libraries
library(afex); library(brms); library(dplyr)

# 80 randomly chosen participants from Study 1; see ? afex::stroop ------
# drawn with set.seed(1); sample(unique(pno), size = 80) from Study 1, and stored here
ids <- c("s1_103", "s1_104", "s1_105", "s1_106", "s1_110", "s1_111", "s1_115", "s1_119",
         "s1_121", "s1_122", "s1_126", "s1_129", "s1_13", "s1_130", "s1_134", "s1_14",
         "s1_146", "s1_148", "s1_156", "s1_160", "s1_162", "s1_163", "s1_165", "s1_166",
         "s1_167", "s1_172", "s1_176", "s1_179", "s1_182", "s1_183", "s1_186", "s1_187",
         "s1_198", "s1_20", "s1_205", "s1_207", "s1_209", "s1_21", "s1_210", "s1_213",
         "s1_215", "s1_216", "s1_217", "s1_218", "s1_22", "s1_220", "s1_222", "s1_223",
         "s1_224", "s1_225", "s1_229", "s1_23", "s1_237", "s1_239", "s1_24", "s1_244",
         "s1_246", "s1_249", "s1_25", "s1_33", "s1_34", "s1_37", "s1_39", "s1_40", "s1_42",
         "s1_43", "s1_44", "s1_45", "s1_51", "s1_68", "s1_7", "s1_70", "s1_73", "s1_74",
         "s1_75", "s1_79", "s1_84", "s1_85", "s1_87", "s1_89")

stroop1 <- subset(afex::stroop, subset = study == 1)
# the stored ids are those drawn with the seed
set.seed(1)
stopifnot(setequal(ids, as.character(sample(unique(stroop1$pno), size = 80))))

data_stroop <- subset(stroop1, subset = pno %in% ids) |>
  # remove useless levels id
  dplyr::mutate(id = droplevels(pno)) |>
  # select relevant variables
  dplyr::select(id, condition, congruency, acc, rt) |>
  # previous cong
  group_by(id, condition) |>
  mutate(prev_cong = lag(congruency)) |>
  ungroup() |>
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
  congruency * condition +
  (1 + prev_cong + congruency + congruency*condition || id)

mod_rt <- brm(formula = formula,
              family = lognormal(),
              cores = 3, chains = 3,
              iter = 8000, seed = 1,
              save_pars = save_pars(all = FALSE, group = FALSE),
              data = data_rt)

save(mod_rt, file = "paper/script/example/mod_rt.rda")

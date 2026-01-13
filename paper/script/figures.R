# Code for manuscript figures docx

pkg = c("dplyr","MASS","tidyr","shape","ggplot2",
        "tibble","furrr","tibble","purrr",
        "mvtnorm", "patchwork","ggpubr")

invisible(sapply(pkg, require, character.only = T))


theme_clean <- function(base_size = 20) {
  list(
    theme_light(base_size = base_size),
    theme(
      axis.line  = element_line(colour = "black", linewidth = 1),
      axis.ticks = element_line(colour = "black", linewidth = 1),
      axis.ticks.length = grid::unit(10, "pt"),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
      # facet strips
      strip.background = element_rect(fill = "grey20", colour = "grey20"),
      strip.text = element_text(colour = "white")
    ),
    guides(
      x = guide_axis(cap = "both"),
      y = guide_axis(cap = "both")
    )
  )
}


# Load simulation results and compute relevant metrics

# Compute FWER
sign_pd = function(data, th){
  data |> 
    group_by(eff, m, s, n, nsim) |> 
    reframe(sign = any(value > th))  
}

sign_pval = function(data, th){
  data |> 
    group_by(eff, m, s, n, nsim) |> 
    reframe(sign = any(value < th))  
}

# Simulation 1
load(file = "paper/script/output/sim1.rda")
sim1 = tibble(sim)
sim_res1 = sim1 |> 
  unnest(res) |> 
  group_by(eff,m, s, n) |> 
  mutate(nsim = 1:n()) |> 
  ungroup() |> 
  unnest(res) |> 
  pivot_longer(c("pval","pval_sidak", "pd", "pd_adj")) |>
  group_by(name) |> 
  nest()

sim_res1$fwer = vector(mode = "list", length = nrow(sim_res1))

for(i in 1:nrow(sim_res1)){
  if(sim_res1$name[i] == "pval" | sim_res1$name[i]== "pval_sidak"){
    sim_res1$fwer[[i]] = sign_pval(sim_res1$data[[i]], 0.05)
  } else {
    sim_res1$fwer[[i]] = sign_pd(sim_res1$data[[i]], .975)}
}


# Simulation 2

# Compute FWER
sign_pd = function(data, th){
  data |> 
    group_by(m,q,eff,nsim) |> 
    reframe(sign = any(value > th))  
}

sign_pval = function(data, th){
  data |> 
    group_by(m,q, eff,nsim) |> 
    reframe(sign = any(value < th))  
}

load(file = "paper/script/output/sim2.rda")
sim2 = tibble(sim)

sim_res2 = sim2 |> 
  unnest(res) |> 
  group_by(m,q,eff) |> 
  mutate(nsim = 1:n()) |> 
  ungroup() |> 
  unnest(res) |> 
  pivot_longer(c("pval","pval_sidak", "pd", "pd_adj")) |>
  group_by(name) |> 
  nest()

sim_res2$fwer = vector(mode = "list", length = nrow(sim_res2))

for(i in 1:nrow(sim_res2)){
  if(sim_res2$name[i] == "pval" | sim_res2$name[i]== "pval_sidak"){
    sim_res2$fwer[[i]] = sign_pval(sim_res2$data[[i]], 0.05)
  } else {
    sim_res2$fwer[[i]] = sign_pd(sim_res2$data[[i]], .975)
  }
}


# Figure 1---------------
set.seed(1)

# Example posterior draws
draws = rnorm(40000, mean = 0.35, sd = 0.25)

# Directional posterior probabilities
p_pos = mean(draws > 0)
p_neg = mean(draws < 0)

pd = max(p_pos, p_neg)

# Density for shaded area plot
dens = density(draws, n = 2000)

df = tibble(x = dens$x, y = dens$y) %>%
  mutate(
    side = if_else(x >= 0, "\u03B8 > 0", "\u03B8 < 0"),
    side = factor(side, levels = c("\u03B8 < 0", "\u03B8 > 0"))
  )

label_txt = paste0(
  "Pr(\u03B8 > 0 | y) = ", round(p_pos, 3), "\n",
  "Pr(\u03B8 < 0 | y) = ", round(p_neg, 3), "\n",
  "pd = ", round(pd, 3)
)

fig1 = ggplot(df, aes(x = x, y = y)) +
  geom_area(aes(fill = side), alpha = 0.8) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  annotate(
    "label",
    x = 1,
    y = max(df$y) * 0.8,
    hjust = 0.5,
    vjust = 0.5,
    label = label_txt,
    size = 3.2
  ) +
  labs(
    x = expression(theta),   # axis label in plotmath (also fine to use "\u03B8")
    y = "Posterior density",
    fill = NULL
  ) +
  scale_fill_manual(
    values = c("\u03B8 < 0" = "deepskyblue", "\u03B8 > 0" = "deepskyblue4")
  ) +
  theme(legend.position = "top")+
  theme_clean(base_size = 12)

ggsave("paper/figures/fig1-posterior-pd.png",  fig1,
       width = 5, height = 3, dpi = 300, bg = "white")


# Figure 2---------------

set.seed(1)

# Case 1: Narrow posterior (excludes zero) 
drawsh1 = rnorm(40000, mean = 0.35, sd = 0.16)

p_pos1 = mean(drawsh1 > 0)
p_neg1 = mean(drawsh1 < 0)
pd1 = max(p_pos1, p_neg1)
p_two_sided1 = 2 * (1 - pd1)

ci_level = 0.95
alpha = (1 - ci_level) / 2
ci1 = quantile(drawsh1, probs = c(alpha, 1 - alpha))
ci_low1 = ci1[1]
ci_high1 = ci1[2]

dens1 = density(drawsh1, n = 2000)
df1 = tibble(x = dens1$x, y = dens1$y) |>
  mutate(
    side = if_else(x >= 0, "θ > 0", "θ < 0"),
    side = factor(side, levels = c("θ < 0", "θ > 0"))
  )

label_txt1 = paste0(
  "Pr(θ > 0 | y) = ", round(p_pos1, 3), "\n",
  "Pr(θ < 0 | y) = ", round(p_neg1, 3), "\n",
  "pd = ", round(pd1, 3), "\n",
  "p = ", signif(p_two_sided1, 3), "\n",
  ci_level*100, "% CI = [", round(ci_low1, 3), ", ", round(ci_high1, 3), "]"
)

h1 = ggplot(df1, aes(x = x, y = y)) +
  geom_area(aes(fill = side), alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_segment(
    aes(x = ci_low1, xend = ci_high1, y = 0, yend = 0),
    linewidth = 3, lineend = "round",
    colour = "black"
  ) +
  annotate(
    "label", x = -.5, y = 1.8,
    hjust = 0.5, vjust = 0.5, label = label_txt1, size = 4.5
  ) +
  labs(x = expression(theta), y = "Posterior density", fill = NULL) +
  scale_fill_manual(values = c("θ < 0" = "deepskyblue2", "θ > 0" = "deepskyblue4")) +
  theme(legend.position = "top") + theme_clean(base_size = 12)+
  scale_x_continuous(breaks = c(-.5,0,.5,1), limits = c(-1,1))

# Case 2: Wide posterior (includes zero) 
draws = rnorm(40000, mean = 0.3, sd = 0.16)

p_pos = mean(draws > 0)
p_neg = mean(draws < 0)
pd = max(p_pos, p_neg)
p_two_sided = 2 * (1 - pd)

ci = quantile(draws, probs = c(alpha, 1 - alpha))
ci_low = ci[1]
ci_high = ci[2]

dens = density(draws, n = 2000)
df = tibble(x = dens$x, y = dens$y) |>
  mutate(
    side = if_else(x >= 0, "θ > 0", "θ < 0"),
    side = factor(side, levels = c("θ < 0", "θ > 0"))
  )

label_txt = paste0(
  "Pr(θ > 0 | y) = ", round(p_pos, 3), "\n",
  "Pr(θ < 0 | y) = ", round(p_neg, 3), "\n",
  "pd = ", round(pd, 3), "\n",
  "p = ", signif(p_two_sided, 3), "\n",
  ci_level*100, "% CI = [", round(ci_low, 3), ", ", round(ci_high, 3), "]"
)

h0 = ggplot(df, aes(x = x, y = y)) +
  geom_area(aes(fill = side), alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, y = 0, yend = 0),
    linewidth = 3, lineend = "round",
    colour = "black"
  ) +
  annotate(
    "label", x = -.5, y = 1.8,
    hjust = 0.5, vjust = 0.5, label = label_txt, size = 4.5
  ) +
  labs(x = expression(theta), y = " ", fill = NULL) +
  scale_fill_manual(values = c("θ < 0" = "deepskyblue", "θ > 0" = "deepskyblue4")) +
  theme(legend.position = "top") + theme_clean(base_size = 12)+
  scale_x_continuous(breaks = c(-.5,0,.5,1), limits = c(-1,1))

fig2 = (h1 | h0) + plot_layout(guides = "collect")&
  theme(legend.position = "top") 

ggsave("paper/figures/fig2-pd-ci.png", fig2,
       width = 10, height = 4, dpi = 300, bg = "white")


# Figure 3---------------

# Distribution of pd under the null 
pval = sim_res1[[2]][[1]]
pd = sim_res1[[2]][[3]]

# Extract relevant info for plotting
df_pval = pval |>
  subset(eff == 0 & m == 1)

bw_pval=0.025
brks_pval=seq(0, 1, by = bw_pval)

pval_den = ggplot(df_pval, aes(x = value)) +
  geom_histogram(
    breaks = brks_pval,
    aes(y = after_stat(density)),alpha = .7, 
    fill = "deeppink4", color = "white",
  ) +
  geom_vline(xintercept = 0.05, size = 1.5, alpha = .8)+
  geom_hline(yintercept = 1, size = 1.5, alpha = .8)+
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "p-value", y = " ") +
  theme_clean(base_size = 18)+
  theme(legend.position = "top")


# Extract relevant info for plotting
df_pd = pd |>
  subset(eff == 0 & m == 1) |>
  mutate(
    n = factor(n, levels = unique(n), labels = paste0("n = ", unique(n))),
    s = factor(s, levels = unique(s), labels = paste0("\u03C4 = ", unique(s)))
  )

bw_pd=0.025
brks_pd=seq(0.5, 1, by = bw_pd)

pd_den = ggplot(df_pd, aes(x = value)) +
  geom_histogram(
    breaks = brks_pd,
    aes(y = after_stat(density)),alpha = .7,
    fill = "deepskyblue4", color = "white",
  ) +
  facet_grid(n~s)+
  geom_vline(xintercept = 0.975,  size = 1.5, alpha = .8)+
  geom_hline(yintercept = 2,  size = 1.5, alpha = .8)+
  coord_cartesian(xlim = c(0.5, 1)) +
  labs(x = "pd", y = " ") +
  theme_clean(base_size = 18)+
  theme(legend.position = "top")

fig3 = (pval_den / pd_den) + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1, 3.5)) 


ggsave("paper/figures/fig3-uniform-pd.png", fig3,
       width = 14, height = 14, dpi = 300, bg = "white")

# Figure 4---------------

# Plot family wise error rate
fig4 = sim_res1 |> 
  unnest(fwer) |> 
  subset(eff == 0 & (name == "pd" | name == "pval"))|>
  group_by(name, m, s,n) |> 
  summarise(FWER = mean(sign)) |> 
  mutate(
    n = factor(n, levels = unique(n), labels = paste0("n = ", unique(n))),
    s = factor(s, levels = unique(s), labels = paste0("\u03C4 = ", unique(s)))
  )|>
  ggplot(aes(x = m, y = FWER, color =name)) +
  geom_line(aes(colour = name), size = 1, show.legend = TRUE) +
  geom_point(aes(colour = name),
             size = 2.5, show.legend = TRUE)+
  scale_colour_manual(values = c(pd = "deepskyblue3", pval = "deeppink3")) +
  facet_grid( n ~ s) +
  geom_hline(yintercept = 0.05, linetype = "dashed")+
  scale_y_continuous(breaks = seq(0,1,.1))+
  scale_x_continuous(breaks = c(1,2,5,10,20))+
  xlab("m")+ labs(colour = "index")+
  theme_clean(base_size = 18)+
  theme(legend.position = "top")

ggsave("paper/figures/fig4-fwer-pd.png", fig4,
       width = 14, height = 10, dpi = 300, bg = "white")

# Figure 5---------------

m = 1:20
sim_k = data.frame(k = sapply(m, function(x) 1 - 0.5^(1/x)), m = m)
fig5 = ggplot(sim_k, aes(y = k, x = m) )+
  geom_line(size = 0.8,color = "deepskyblue4")+
  geom_point(size = 2, color = "black", alpha = .6)+
  scale_y_continuous(breaks = seq(0,0.5,0.1))+
  scale_x_continuous(breaks = round(seq(1,20,1)))+
  theme_clean(base_size = 14)+
  theme(legend.position = "top")

ggsave("paper/figures/fig5-prior-adj.png", fig5,
       width = 7, height = 4, dpi = 300, bg = "white")

# Figure 6---------------

# Plot family wise error rate
fig6 = sim_res1 |> 
  unnest(fwer) |> 
  subset(name == "pd_adj" | name == "pval_sidak")|>
  group_by(name, m, s,n) |> 
  summarise(FWER = mean(sign)) |> 
  mutate(
    n = factor(n, levels = unique(n), labels = paste0("n = ", unique(n))),
    s = factor(s, levels = unique(s), labels = paste0("\u03C4 = ", unique(s)))
  )|>
  ggplot(aes(x = m, y = FWER, color =name)) +
  geom_line(aes(colour = name), size = 0.8, show.legend = TRUE) +
  geom_point(aes(colour = name),
             size = 2, show.legend = TRUE)+
  scale_colour_manual(values = c(pd_adj = "deepskyblue3",
                                 pval_sidak = "deeppink3")) +
  facet_grid( n ~ s) +
  geom_hline(yintercept = 0.05, linetype = "dashed")+
  scale_y_continuous(breaks = seq(0,.05,.01))+
  scale_x_continuous(breaks = c(1,2,5,10,20))+
  xlab("m")+labs(colour = "index")+
  theme_clean(base_size = 14)+
  theme(legend.position = "top")

ggsave("paper/figures/fig6-fwer-adj.png", fig6,
       width = 10, height = 8, dpi = 300, bg = "white")

# Figure 8---------------

# Plot family wise error rate
pnull = sim_res2 |> 
  unnest(fwer) |> 
  subset(eff == 0 & (name == "pd_adj" | name == "pval_sidak"))|>
  group_by(name, m, q) |> 
  summarise(FWER = mean(sign)) |>
  ggplot(aes(x = m, y = FWER, color =name)) +
  geom_line(aes(colour = name), size = 1, show.legend = TRUE) +
  geom_point(aes(colour = name),
             size = 2.5, show.legend = TRUE)+
  scale_colour_manual(values = c(pd_adj = "deepskyblue3",
                                 pval_sidak = "deeppink3")) +
  facet_grid(~q,
             labeller = as_labeller(c(
               `0.1` = "Pr(all null) = 0.1",
               `0.2` = "Pr(all null) = 0.2",
               `0.3` = "Pr(all null) = 0.3",
               `0.4` = "Pr(all null) = 0.4",
               `0.5` = "Pr(all null) = 0.5"))) +
  scale_y_continuous(breaks = seq(0,.25,.05))+
  scale_x_continuous(breaks = c(1,2,3,4,5,10))+
  xlab("m")+
  geom_hline(yintercept = 0.05, linetype = "dashed")+
  labs(colour = "index")+
  theme_clean(base_size = 14)+ theme(legend.position = "top")

ppower = sim_res2 |> 
  unnest(fwer) |> 
  subset(eff == 0.2 & (name == "pd_adj" | name == "pval_sidak"))|>
  group_by(name, q,m) |> 
  summarise(POWER = mean(sign)) |>
  ggplot(aes(x = m, y = POWER, color =name)) +
  geom_line(aes(colour = name), size = 1, show.legend = TRUE) +
  geom_point(aes(colour = name),
             size = 2.5, show.legend = TRUE)+
  scale_colour_manual(values = c(pd_adj = "deepskyblue3",
                                 pval_sidak = "deeppink3")) +
  facet_grid(~q,
             labeller = as_labeller(c(
               `0.1` = "Pr(all null) = 0.1",
               `0.2` = "Pr(all null) = 0.2",
               `0.3` = "Pr(all null) = 0.3",
               `0.4` = "Pr(all null) = 0.4",
               `0.5` = "Pr(all null) = 0.5"))) +
  #scale_y_continuous(breaks = seq(0,.25,.025))+
  scale_x_continuous(breaks = c(1,2,3,4,5,10))+
  xlab("m")+
  labs(colour = "index")+
  theme_clean(base_size = 14)+ theme(legend.position = "top")

diff= sim_res2 |> 
  unnest(fwer) |> 
  subset(name == "pd_adj" | name == "pval_sidak")|>
  group_by(name, q,m,eff) |> 
  summarise(fwer = mean(sign)) |>
  ungroup()|>
  group_by(q,m)|>
  summarise(diff_fwer = fwer[name=="pd_adj"& eff ==0]-
              fwer[name=="pval_sidak"& eff ==0],
            diff_pw = fwer[name=="pd_adj"& eff ==0.2]-
              fwer[name=="pval_sidak"& eff ==0.2],
            diff_diff = diff_pw-diff_fwer) |>
  ggplot(aes(x = m, y = diff_fwer)) +
  geom_line(aes(y = diff_fwer, colour = "ΔFWER"), linewidth = 1) +
  geom_point(aes(y = diff_fwer, colour = "ΔFWER"), size = 2.5) +
  geom_line(aes(y = diff_pw, colour = "ΔPower"), linewidth = 1) +
  geom_point(aes(y = diff_pw, colour = "ΔPower"), size = 2.5) +
  scale_colour_manual(
    name = NULL,
    values = c("ΔFWER" = "red4",
               "ΔPower" = "green4")
  ) +
  facet_grid(~q,
             labeller = as_labeller(c(
               `0.1` = "Pr(all null) = 0.1",
               `0.2` = "Pr(all null) = 0.2",
               `0.3` = "Pr(all null) = 0.3",
               `0.4` = "Pr(all null) = 0.4",
               `0.5` = "Pr(all null) = 0.5"))) +
  scale_y_continuous(breaks = round(seq(-0.1,.2,.1),2), limits = c(-.15,.25))+
  scale_x_continuous(breaks = c(1,2,3,4,5,10))+
  xlab("m")+ylab("Δ = pd - p-value")+
  labs(colour = "index")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_clean(base_size = 14)+ theme(legend.position = "top")

fig7 = (pnull/ppower/diff) + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(.8, 1,1.5)) 


ggsave("paper/figures/fig7-fwer-power.png", fig7,
       width = 10, height = 12, dpi = 300, bg = "white")





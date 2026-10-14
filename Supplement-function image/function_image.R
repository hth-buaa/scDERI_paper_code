###########################
#                         #
#      Function Plot      #
#                         #
###########################



setwd("/mnt/i/scDEDM_Paper/Supplement-function image/")
library(ggplot2)



########## 1. Hill Function H Image ##########
x = seq(0, 7, length.out = 300)
hill_func = function(x, K, n) x^n / (K + x^n)

df_hill = data.frame(
  x = rep(x, 4),
  y = c(hill_func(x, 2, 2), hill_func(x, 2, 1), hill_func(x, 2, 0.5), rep(1, length(x))),
  func = rep(c("H1", "H2", "H3", "asymptote"), each = length(x))
)
df_hill$func = factor(df_hill$func, levels = c("H1", "H2", "H3", "asymptote"))

color_hill = c("H1" = "darkblue", "H2" = "darkgreen", "H3" = "darkorange", "asymptote" = "black")
linetype_hill = c("H1" = "solid", "H2" = "solid", "H3" = "solid", "asymptote" = "dashed")
legend_hill = c(
  expression(K[H] == 2 * "," ~ n == 2),
  expression(K[H] == 2 * "," ~ n == 1),
  expression(K[H] == 2 * "," ~ n == 0.5),
  expression(y == 1)
)

plot_H = ggplot(df_hill, aes(x, y, color = func, linetype = func)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = color_hill, labels = legend_hill) +
  scale_linetype_manual(values = linetype_hill, labels = legend_hill) +
  labs(x = "x", y = "y", title = "Hill Function") +
  coord_cartesian(xlim = c(0, 7)) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "right", legend.title = element_blank(),
    legend.text = element_text(size = 9), legend.key.width = unit(1, "cm")
  )
print(plot_H)
ggsave(filename = "./Hill_function.png", plot = plot_H, width = 5, height = 3, dpi = 300)


########## 2. Derivation of Impact Scale S ##########
##### 2.1 h0 Correlation Function Image #####
x = seq(0, 3, length.out = 300)
df_h0 = data.frame(
  x = rep(x, 6),
  y = c(x, rep(1, length(x)), exp(-x), -exp(-x), x*exp(-x), (1-x)*exp(-x)),
  func = rep(c("h01", "h01_prime", "h02", "h02_prime", "h0", "h0_prime"), each = length(x))
)
df_h0$func = factor(df_h0$func, levels = c("h01", "h01_prime", "h02", "h02_prime", "h0", "h0_prime"))

color_h0 = c(
  "h01" = "darkblue", "h01_prime" = "darkblue",
  "h02" = "darkgreen", "h02_prime" = "darkgreen",
  "h0" = "darkorange", "h0_prime" = "darkorange"
)
linetype_h0 = c(
  "h01" = "solid", "h01_prime" = "dashed",
  "h02" = "solid", "h02_prime" = "dashed",
  "h0" = "solid", "h0_prime" = "dashed"
)
legend_h0 = c(
  expression(h[0*","*1](x)), expression(h[0*","*1]^{"'"}*(x)),
  expression(h[0*","*2](x)), expression(h[0*","*2]^{"'"}*(x)),
  expression(h[0](x)), expression(h[0]^{"'"}*(x))
)

plot_h0 = ggplot(df_h0, aes(x, y, color = func, linetype = func)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = color_h0, labels = legend_h0) +
  scale_linetype_manual(values = linetype_h0, labels = legend_h0) +
  labs(x = "x", y = "y", title = expression(h[0]~" Correlation Function")) +
  coord_cartesian(xlim = c(0, 3)) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "right", legend.title = element_blank(),
    legend.text = element_text(size = 9), legend.key.width = unit(1, "cm")
  )
print(plot_h0)
ggsave(filename = "./h0.png", plot = plot_h0, width = 5, height = 3, dpi = 300)


##### 2.2 h1 Correlation Function Image #####
x = seq(0, 1, length.out = 300)
x0 = 0.6
df_h1 = data.frame(
  x = rep(x, 3),
  y = c(
    (x/x0)*exp(1 - x/x0),
    ((1-x)/(1-x0))*exp(1 - (1-x)/(1-x0)),
    (x*(1-x)/(x0*(1-x0)))*exp(2 - x/x0 - (1-x)/(1-x0))
  ),
  func = rep(c("h11", "h12", "h1"), each = length(x))
)
df_h1$func = factor(df_h1$func, levels = c("h11", "h12", "h1"))

color_h1 = c("h11" = "darkblue", "h12" = "darkgreen", "h1" = "darkorange")
linetype_h1 = c("h11" = "solid", "h12" = "solid", "h1" = "solid")
legend_h1 = c(expression(h[1*","*1](x)), expression(h[1*","*2](x)), expression(h[1](x)))

plot_h1 = ggplot(df_h1, aes(x, y, color = func, linetype = func)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = color_h1, labels = legend_h1) +
  scale_linetype_manual(values = linetype_h1, labels = legend_h1) +
  labs(x = "x", y = "y", title = expression(h[1]~" Correlation Function"~(x[0]==0.6))) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1, "cm")
  )
print(plot_h1)
ggsave(filename = "./h1.png", plot = plot_h1, width = 5, height = 3, dpi = 300)


##### 2.3 h2 Correlation Function Image #####
x = seq(0, 1, length.out = 300)
h2 = function(x, x0, n) ((x*(1-x)/(x0*(1-x0)))*exp(2 - x/x0 - (1-x)/(1-x0)))^n
h2_prime = function(x, x0, n) {
  numerator = n * x^(n-1) * (1-x)^(n-1) * (x0 - x) * ((1-x0)*(1-x) + x0*x)
  denominator = x0^(n+1) * (1-x0)^(n+1) * exp(n*(2*x0 - 1)*(x0 - x)/(x0*(1-x0)))
  numerator / denominator
}

params = list(
  list(x0 = 0.6, n = 0.5, color = "darkblue"),
  list(x0 = 0.75, n = 0.8, color = "darkgreen"),
  list(x0 = 0.25, n = 0.8, color = "darkorange"),
  list(x0 = 0.6, n = 1, color = "darkred"),
  list(x0 = 0.6, n = 2, color = "darkgrey")
)

df_h2 = data.frame()
for (i in seq_along(params)) {
  p = params[[i]]
  y_h2 = h2(x, p$x0, p$n)
  y_h2p = h2_prime(x, p$x0, p$n)
  y_h2p[y_h2p > 3 | y_h2p < -4] = NA

  df_h2 = rbind(
    df_h2,
    data.frame(x = x, y = y_h2, func = paste0("h2_", i), type = "func"),
    data.frame(x = x, y = y_h2p, func = paste0("h2p_", i), type = "prime")
  )
}
df_h2$func = factor(df_h2$func, levels = c(paste0("h2_", 1:5), paste0("h2p_", 1:5)))

color_h2 = setNames(
  rep(sapply(params, function(p) p$color), 2),
  c(paste0("h2_", 1:5), paste0("h2p_", 1:5))
)
linetype_h2 = setNames(
  rep(c("solid", "dashed"), each = 5),
  c(paste0("h2_", 1:5), paste0("h2p_", 1:5))
)
legend_h2 = c(
  expression(h[2](x) * "," ~ x[0] == 0.6 * "," ~ n == 0.5),
  expression(h[2](x) * "," ~ x[0] == 0.75 * "," ~ n == 0.8),
  expression(h[2](x) * "," ~ x[0] == 0.25 * "," ~ n == 0.8),
  expression(h[2](x) * "," ~ x[0] == 0.6 * "," ~ n == 1),
  expression(h[2](x) * "," ~ x[0] == 0.6 * "," ~ n == 2),
  expression(h[2]^{"'"}*(x) * "," ~ x[0] == 0.6 * "," ~ n == 0.5),
  expression(h[2]^{"'"}*(x) * "," ~ x[0] == 0.75 * "," ~ n == 0.8),
  expression(h[2]^{"'"}*(x) * "," ~ x[0] == 0.25 * "," ~ n == 0.8),
  expression(h[2]^{"'"}*(x) * "," ~ x[0] == 0.6 * "," ~ n == 1),
  expression(h[2]^{"'"}*(x) * "," ~ x[0] == 0.6 * "," ~ n == 2)
)

plot_h2 = ggplot(df_h2, aes(x, y, color = func, linetype = func)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  scale_color_manual(values = color_h2, labels = legend_h2) +
  scale_linetype_manual(values = linetype_h2, labels = legend_h2) +
  labs(x = "x", y = "y", title = expression(h[2]~" Correlation Function")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-4, 3)) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
    legend.position = "right", legend.title = element_blank(),
    legend.text = element_text(size = 9), legend.key.width = unit(1, "cm")
  )
print(plot_h2)
ggsave(filename = "./h2.png", plot = plot_h2, width = 5, height = 3, dpi = 300)



########## 3. Impact Efficiency R ##########
theta = seq(0, 1, length.out = 300)
tau = 0.4
r1 = 200; r2 = 1; r3 = 0.5; r4 = 1000; r5 = 1

term_norm_1 = 1 / (r1 * tau)
term_norm_2 = r5 / ((1 + 10^(r2 * (r3 - 1))) * (1 + 10^(r4 * r2 * (tau - 1))))
r6 = term_norm_1 + term_norm_2

R_func = function(theta) {
  term1 = theta / (r1 * tau)
  denominator_part1 = 1 + 10^(r2 * (r3 - theta))
  denominator_part2 = 1 + 10^(r4 * r2 * (tau - theta))
  term2 = r5 / (denominator_part1 * denominator_part2)
  (term1 + term2) / r6
}

df_R = data.frame(x = theta, y = R_func(theta))

color_R = c("R" = "darkblue", "vline" = "black")
linetype_R = c("R" = "solid", "vline" = "dashed")
legend_R = c(expression(R(θ)), expression(x == τ))

plot_R = ggplot(df_R, aes(x, y)) +
  geom_line(aes(color = "R", linetype = "R"), linewidth = 0.7) +
  geom_vline(aes(xintercept = tau, color = "vline", linetype = "vline")) +
  scale_color_manual(values = color_R, labels = legend_R, breaks = c("R", "vline"), name = "") +
  scale_linetype_manual(values = linetype_R, labels = legend_R, breaks = c("R", "vline"), name = "") +
  guides(color = guide_legend(), linetype = guide_legend()) +
  labs(x = expression(θ), y = expression(R),
       title = expression("Impact Efficiency R " ~ (τ==0.4 ~ "," ~ r[1]==200 ~ "," ~ r[2]==1 ~ "," ~ r[3]==0.5 ~ "," ~ r[4]==1000 ~ "," ~ r[5]==1))) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
    legend.position = "right", legend.title = element_blank(),
    legend.text = element_text(size = 9), legend.key.width = unit(1, "cm")
  )
print(plot_R)
ggsave(filename = "./Impact_Efficiency_R.png", plot = plot_R, width = 5, height = 3, dpi = 300)

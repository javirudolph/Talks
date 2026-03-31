# ============================================================
#  Individual variation in movement & disease transmission
#  Two animals, same mean step length, different tail shape
#  Step lengths ~ Gamma(shape, rate)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(42)

# ── Parameters ───────────────────────────────────────────────

mean_step  <- 5       # same mean for both animals (e.g. km)
shape_avg  <- 10       # average mover: concentrated, light tail
shape_tail <- 0.5     # heavy-tailed mover: dispersed, fat tail

# For Gamma: mean = shape/rate  =>  rate = shape/mean
rate_avg   <- shape_avg  / mean_step
rate_tail  <- shape_tail / mean_step

n_steps    <- 500     # steps per trajectory
d_star     <- 15      # transmission threshold (landscape boundary)

# ── 1. Draw step lengths ─────────────────────────────────────

steps_avg  <- rgamma(n_steps, shape = shape_avg,  rate = rate_avg)
steps_tail <- rgamma(n_steps, shape = shape_tail, rate = rate_tail)

cat("Average mover   — mean:", round(mean(steps_avg),  2),
    "  max:", round(max(steps_avg),  2), "\n")
cat("Heavy-tail mover — mean:", round(mean(steps_tail), 2),
    "  max:", round(max(steps_tail), 2), "\n")

# ── 2. Random walk trajectories ──────────────────────────────

make_walk <- function(steps) {
  angles <- runif(length(steps), 0, 2 * pi)
  x <- cumsum(steps * cos(angles))
  y <- cumsum(steps * sin(angles))
  data.frame(x = c(0, x), y = c(0, y),
             step = 0:length(steps))
}

walk_avg  <- make_walk(steps_avg)  %>% mutate(animal = "Average mover")
walk_tail <- make_walk(steps_tail) %>% mutate(animal = "Heavy-tail mover")
walks     <- bind_rows(walk_avg, walk_tail)

# ── 3. Step length distributions ─────────────────────────────

x_grid <- seq(0, 40, length.out = 500)

dens_df <- bind_rows(
  data.frame(x     = x_grid,
             dens  = dgamma(x_grid, shape = shape_avg,  rate = rate_avg),
             animal = "Average mover"),
  data.frame(x     = x_grid,
             dens  = dgamma(x_grid, shape = shape_tail, rate = rate_tail),
             animal = "Heavy-tail mover")
)

# ── 4. Transmission probability P(step >= d*) ────────────────

# Analytically: P(X >= d*) = 1 - pgamma(d*, shape, rate)
p_avg  <- 1 - pgamma(d_star, shape = shape_avg,  rate = rate_avg)
p_tail <- 1 - pgamma(d_star, shape = shape_tail, rate = rate_tail)

cat("\nP(step >= d* =", d_star, "):\n")
cat("  Average mover:    ", round(p_avg,  5), "\n")
cat("  Heavy-tail mover: ", round(p_tail, 4), "\n")

# ── 5. Transmission events over infectious period ────────────
#
# For each simulated step: did the animal cross d*?
# Over N_infectious steps, count how many crossings occur.
# Repeat across many individuals drawn from each population.

N_infectious <- 20    # steps taken while infectious
N_individuals <- 2000 # individuals in population

sim_transmission <- function(shape, rate, n_ind, n_inf, threshold) {
  vapply(seq_len(n_ind), function(i) {
    steps <- rgamma(n_inf, shape = shape, rate = rate)
    sum(steps >= threshold)   # number of transmission events
  }, integer(1L))
}

trans_avg  <- sim_transmission(shape_avg,  rate_avg,  N_individuals,
                               N_infectious, d_star)
trans_tail <- sim_transmission(shape_tail, rate_tail, N_individuals,
                               N_infectious, d_star)

trans_df <- bind_rows(
  data.frame(events = trans_avg,  animal = "Average mover"),
  data.frame(events = trans_tail, animal = "Heavy-tail mover")
)

cat("\nTransmission events over", N_infectious, "infectious steps:\n")
cat("  Average mover    — mean:", round(mean(trans_avg),  3),
    "  P(>=1 event):", round(mean(trans_avg  >= 1), 3), "\n")
cat("  Heavy-tail mover — mean:", round(mean(trans_tail), 3),
    "  P(>=1 event):", round(mean(trans_tail >= 1), 3), "\n")

# ── 6. EVT: fit GPD to threshold exceedances ─────────────────

library(evd)

all_steps <- c(steps_avg, steps_tail)
u         <- quantile(all_steps, 0.90)   # 90th percentile as threshold
exceedances <- all_steps[all_steps > u] - u

fit   <- fpot(all_steps, threshold = u, std.err = FALSE)
xi    <- round(fit$estimate["shape"], 3)
scale <- round(fit$estimate["scale"], 3)

cat("\nGPD fit to pooled step lengths (threshold =", round(u, 2), "):\n")
cat("  xi (shape) =", xi, "\n")
cat("  scale      =", scale, "\n")

# ── 7. Plots ──────────────────────────────────────────────────

pal <- c("Average mover"    = "#2B7BB9",
         "Heavy-tail mover" = "#C1440E")

# Panel A: step length distributions
pA <- ggplot(dens_df, aes(x = x, y = dens,
                           colour = animal, fill = animal)) +
  geom_area(alpha = 0.25, position = "identity") +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = d_star, linetype = "dashed",
             colour = "grey30", linewidth = 0.8) +
  annotate("text", x = d_star + 0.5, y = Inf,
           label = paste0("d* = ", d_star),
           hjust = 0, vjust = 1.5, size = 3.5, colour = "grey30") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  coord_cartesian(xlim = c(0, 35)) +
  labs(title    = "Step length distributions",
       subtitle = paste0("Same mean (", mean_step,
                         ") — different tail shape"),
       x = "Step length", y = "Density",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

# Panel B: random walk trajectories
pB <- ggplot(walks, aes(x = x, y = y,
                         colour = animal, group = animal)) +
  geom_path(alpha = 0.6, linewidth = 0.5) +
  geom_point(data = walks %>% filter(step == 0),
             size = 3, shape = 21,
             aes(fill = animal), colour = "white") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  coord_equal() +
  labs(title    = "Movement trajectories",
       subtitle = paste0(n_steps, " steps each"),
       x = "x", y = "y",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

# Panel C: transmission events distribution
pC <- ggplot(trans_df, aes(x = events, fill = animal,
                            colour = animal)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, alpha = 0.55,
                 position = "identity") +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(title    = "Transmission events per individual",
       subtitle = paste0("Over ", N_infectious,
                         " infectious steps  |  d* = ", d_star),
       x = "Number of transmission events",
       y = "Density",
       fill = NULL, colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

# Panel D: exceedances + GPD fit
exc_df <- data.frame(x = exceedances)
x_exc  <- seq(0, max(exceedances), length.out = 300)
gpd_df <- data.frame(
  x    = x_exc,
  dens = dgpd(x_exc, loc = 0, scale = scale, shape = xi)
)

pD <- ggplot(exc_df, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 fill = adjustcolor("#7F77DD", 0.5),
                 colour = "white") +
  geom_line(data = gpd_df, aes(x = x, y = dens),
            colour = "grey20", linewidth = 1.0) +
  annotate("text", x = max(exceedances) * 0.6,
           y = max(dgpd(x_exc, 0, scale, xi)) * 0.8,
           label = paste0("GPD fit\n\u03be = ", xi),
           size = 3.8, colour = "grey20") +
  labs(title    = "Threshold exceedances (pooled steps)",
       subtitle = paste0("POT threshold u = ", round(u, 1),
                         "  (90th percentile)"),
       x = "Excess step length  (step \u2212 u)",
       y = "Density") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

# ── Assemble ──────────────────────────────────────────────────

(pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Individual variation in movement drives transmission heterogeneity",
    subtitle = paste0(
      "Same mean step length (", mean_step, ") — ",
      "average mover: Gamma(shape=", shape_avg, ")  |  ",
      "heavy-tail mover: Gamma(shape=", shape_tail, ")"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 10)
    )
  )

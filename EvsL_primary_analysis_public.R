library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(12345)

# Parameters

t_max <- 11
S_min <- 0
S_max <- 9
mu0 <- 4.5
sigmaS <- 2.0
beta_S <- 3.5

alpha0 <- -4.0
alphaL <- 0.5
alpha_t <- 0.05
gamma_delay <- 0.8
delta_benefit <- 0.3

# Data-generating mechanism

gen_L0 <- function(n) {
  rnorm(n, mean = 0, sd = 1)
}

generate_S_scenario_A <- function(n, S_min, S_max, mu0, sigmaS) {
  times_S <- S_min:S_max
  w_base <- exp(-(times_S - mu0)^2 / (2 * sigmaS^2))
  w_base <- w_base / sum(w_base)
  sample(times_S, size = n, replace = TRUE, prob = w_base)
}

generate_S_scenario_B <- function(L0, S_min, S_max, mu0, sigmaS, beta_S) {
  times_S <- S_min:S_max
  vapply(L0, function(x) {
    mu_i <- mu0 - beta_S * x
    w_i <- exp(-(times_S - mu_i)^2 / (2 * sigmaS^2))
    w_i <- w_i / sum(w_i)
    sample(times_S, size = 1, replace = TRUE, prob = w_i)
  }, integer(1))
}

logit_to_p <- function(x) {
  exp(x) / (1 + exp(x))
}

simulate_event_time <- function(L0_i, S_i,
                                sim_type = 1,
                                t_max = 11,
                                alpha0 = -4.0,
                                alphaL = 0.5,
                                alpha_t = 0.05,
                                gamma_delay = 0.8,
                                delta_benefit = 0.3) {
  event_time <- NA_integer_

  for (t in 0:t_max) {
    lp <- alpha0 + alphaL * L0_i + alpha_t * t

    if (sim_type == 0) {
      if (t >= S_i) {
        time_on_treat <- t - S_i + 1
        lp <- lp - delta_benefit * time_on_treat
      }
    } else if (sim_type == 1) {
      if (t < S_i) {
        lp <- lp + gamma_delay
      }
    } else if (sim_type == 2) {
      if (t < S_i) {
        lp <- lp + gamma_delay
      } else {
        time_on_treat <- t - S_i + 1
        lp <- lp - delta_benefit * time_on_treat
      }
    } else if (sim_type == 3) {
      lp <- lp
    } else {
      stop("sim_type must be 0, 1, 2, or 3.")
    }

    p_t <- logit_to_p(lp)
    if (runif(1) < p_t) {
      event_time <- t
      break
    }
  }

  if (is.na(event_time)) {
    event_time <- t_max + 1
  }

  event_time
}

simulate_dataset <- function(n,
                             sim_type = 1,
                             scenario = c("A", "B"),
                             t_max = 11,
                             S_min = 0,
                             S_max = 9) {
  scenario <- match.arg(scenario)
  L0 <- gen_L0(n)

  if (scenario == "A") {
    S <- generate_S_scenario_A(
      n = n,
      S_min = S_min,
      S_max = S_max,
      mu0 = mu0,
      sigmaS = sigmaS
    )
  } else {
    S <- generate_S_scenario_B(
      L0 = L0,
      S_min = S_min,
      S_max = S_max,
      mu0 = mu0,
      sigmaS = sigmaS,
      beta_S = beta_S
    )
  }

  event_time <- integer(n)
  status <- integer(n)

  for (i in seq_len(n)) {
    et <- simulate_event_time(
      L0_i = L0[i],
      S_i = S[i],
      sim_type = sim_type,
      t_max = t_max,
      alpha0 = alpha0,
      alphaL = alphaL,
      alpha_t = alpha_t,
      gamma_delay = gamma_delay,
      delta_benefit = delta_benefit
    )
    event_time[i] <- et
    status[i] <- as.integer(et <= t_max)
  }

  max_treat_time <- t_max - S + 1
  observed_treat_time <- pmax(0, pmin(event_time, t_max) - S + 1)
  observed_treat_time[event_time <= S] <- 0

  data.frame(
    id = 1:n,
    L0 = L0,
    S = S,
    event_time = event_time,
    status = status,
    sim_type = sim_type,
    scenario = scenario,
    max_treat_time = max_treat_time,
    observed_treat_time = observed_treat_time
  )
}

# Overlap of baseline covariate distributions

calc_ovl_L0 <- function(dat, e_start, e_end, l_start, l_end) {
  L0_early <- dat$L0[dat$S >= e_start & dat$S <= e_end]
  L0_late <- dat$L0[dat$S >= l_start & dat$S <= l_end]

  if (length(L0_early) < 10 || length(L0_late) < 10) {
    warning("Too few observations in the Early or Late group; OVL set to NA.")
    return(NA_real_)
  }

  rng <- range(c(L0_early, L0_late))
  densE <- stats::density(L0_early, from = rng[1], to = rng[2], n = 512)
  densL <- stats::density(L0_late, from = rng[1], to = rng[2], n = 512)
  dx <- densE$x[2] - densE$x[1]
  sum(pmin(densE$y, densL$y)) * dx
}

# Clone-censor-weight implementation

make_clones_with_censor <- function(dat,
                                    e_start, e_end,
                                    l_start, l_end,
                                    t_max = 11) {
  dat_early <- dat
  dat_early$arm <- "early"

  dat_late <- dat
  dat_late$arm <- "late"

  clone_dat <- rbind(dat_early, dat_late)
  clone_dat$id_clone <- paste0(clone_dat$id, "_", clone_dat$arm)

  calc_censor_time_one <- function(S_i, arm,
                                   e_start, e_end,
                                   l_start, l_end,
                                   t_max) {
    censor_time <- t_max + 1

    if (arm == "early") {
      a <- e_start
      b <- e_end
      if (S_i < a) {
        censor_time <- S_i
      } else if (S_i > b) {
        censor_time <- b + 1
      }
    } else if (arm == "late") {
      c <- l_start
      d <- l_end
      if (S_i < c) {
        censor_time <- S_i
      } else if (S_i > d) {
        censor_time <- d + 1
      }
    }

    censor_time
  }

  clone_dat$censor_time <- mapply(
    FUN = calc_censor_time_one,
    S_i = clone_dat$S,
    arm = clone_dat$arm,
    MoreArgs = list(
      e_start = e_start,
      e_end = e_end,
      l_start = l_start,
      l_end = l_end,
      t_max = t_max
    )
  )

  clone_dat$t_star <- pmin(clone_dat$event_time, clone_dat$censor_time, t_max + 1)
  clone_dat$event_ccw <- as.integer(
    clone_dat$event_time <= t_max & clone_dat$event_time <= clone_dat$censor_time
  )

  clone_dat
}

estimate_ipcw_weights_stab <- function(clone_dat, t_max = 11, cap_w = 20) {
  pp <- lapply(seq_len(nrow(clone_dat)), function(i) {
    r <- clone_dat[i, ]
    ts <- 0:(r$t_star - 1)
    if (length(ts) == 0) return(NULL)

    data.frame(
      id_clone = r$id_clone,
      arm = r$arm,
      L0 = r$L0,
      t = ts,
      cens_this = as.integer(ts == (r$censor_time - 1) & r$censor_time <= t_max)
    )
  }) |>
    bind_rows() |>
    mutate(t_sq = t^2)

  fit_E <- glm(
    cens_this ~ L0 + t + t_sq,
    data = filter(pp, arm == "early"),
    family = binomial()
  )

  fit_L <- glm(
    cens_this ~ L0 + t + t_sq,
    data = filter(pp, arm == "late"),
    family = binomial()
  )

  pp$pcens_hat <- ifelse(
    pp$arm == "early",
    predict(fit_E, pp, type = "response"),
    predict(fit_L, pp, type = "response")
  )

  pp_means <- pp |>
    group_by(arm, t) |>
    summarise(pcens_mean = mean(pcens_hat), .groups = "drop")

  pp <- left_join(pp, pp_means, by = c("arm", "t"))

  G_df <- pp |>
    group_by(id_clone) |>
    summarise(
      G_hat = prod(1 - pcens_hat),
      G_num = prod(1 - pcens_mean),
      .groups = "drop"
    )

  clone_dat |>
    left_join(G_df, by = "id_clone") |>
    mutate(w_ipcw = pmin(G_num / G_hat, cap_w))
}

summarize_risk_rr_ess <- function(clone_weighted) {
  risk_tab <- clone_weighted |>
    group_by(arm) |>
    summarise(
      risk_hat = sum(w_ipcw * event_ccw) / sum(w_ipcw),
      ESS = (sum(w_ipcw)^2) / sum(w_ipcw^2),
      n_raw = n(),
      w_mean = mean(w_ipcw),
      w_max = max(w_ipcw),
      .groups = "drop"
    )

  risk_E <- risk_tab$risk_hat[risk_tab$arm == "early"]
  risk_L <- risk_tab$risk_hat[risk_tab$arm == "late"]
  RR <- as.numeric(risk_E / risk_L)

  list(risk_table = risk_tab, RR = RR)
}

analyze_one_window <- function(dat,
                               e_start, e_end,
                               l_start, l_end,
                               t_max = 11) {
  ovl_L0 <- calc_ovl_L0(dat, e_start, e_end, l_start, l_end)

  clone_dat <- make_clones_with_censor(
    dat,
    e_start = e_start,
    e_end = e_end,
    l_start = l_start,
    l_end = l_end,
    t_max = t_max
  )

  clone_w <- estimate_ipcw_weights_stab(clone_dat, t_max = t_max, cap_w = 20)
  risk_res <- summarize_risk_rr_ess(clone_w)

  list(
    window = c(
      e_start = e_start,
      e_end = e_end,
      l_start = l_start,
      l_end = l_end
    ),
    OVL_L0 = ovl_L0,
    risk_tab = risk_res$risk_table,
    RR = risk_res$RR
  )
}

# Window families

make_windows_sliding <- function() {
  data.frame(
    family = "Sliding",
    window_id = 1:5,
    e_start = 0:4,
    e_end = 0:4 + 2,
    l_start = 0:4 + 3,
    l_end = 0:4 + 5
  )
}

make_windows_late_disjoint <- function() {
  L_end <- 5:9
  data.frame(
    family = "Late_disjoint",
    window_id = seq_along(L_end),
    e_start = 0,
    e_end = 2,
    l_start = 3,
    l_end = L_end
  )
}

make_windows_late_nested <- function() {
  L_end <- 3:9
  data.frame(
    family = "Late_nested",
    window_id = seq_along(L_end),
    e_start = 0,
    e_end = 2,
    l_start = 0,
    l_end = L_end
  )
}

make_windows_gap <- function() {
  L_start <- 3:7
  data.frame(
    family = "Gap_expansion",
    window_id = seq_along(L_start),
    e_start = 0,
    e_end = 2,
    l_start = L_start,
    l_end = L_start + 2
  )
}

make_windows_overlap <- function() {
  data.frame(
    family = "Overlap_crossing",
    window_id = 1:3,
    e_start = c(0, 1, 2),
    e_end = c(4, 5, 6),
    l_start = c(5, 4, 3),
    l_end = c(9, 8, 7)
  )
}

windows_all <- bind_rows(
  make_windows_sliding(),
  make_windows_late_disjoint(),
  make_windows_late_nested(),
  make_windows_gap(),
  make_windows_overlap()
) |>
  mutate(window_name = paste0(family, "_", window_id))

# Generate datasets used in the primary analysis

n_sample <- 5000

sim_datasets <- list(
  Sim0_A = simulate_dataset(n = n_sample, sim_type = 0, scenario = "A", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim0_B = simulate_dataset(n = n_sample, sim_type = 0, scenario = "B", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim1_A = simulate_dataset(n = n_sample, sim_type = 1, scenario = "A", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim1_B = simulate_dataset(n = n_sample, sim_type = 1, scenario = "B", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim2_A = simulate_dataset(n = n_sample, sim_type = 2, scenario = "A", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim2_B = simulate_dataset(n = n_sample, sim_type = 2, scenario = "B", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim3_A = simulate_dataset(n = n_sample, sim_type = 3, scenario = "A", t_max = t_max, S_min = S_min, S_max = S_max),
  Sim3_B = simulate_dataset(n = n_sample, sim_type = 3, scenario = "B", t_max = t_max, S_min = S_min, S_max = S_max)
)

# Run the primary analysis across all window families

all_rr_list <- list()

for (sim_name in names(sim_datasets)) {
  dat <- sim_datasets[[sim_name]]
  message("Running windows for: ", sim_name)

  res_each_sim <- lapply(seq_len(nrow(windows_all)), function(i) {
    w <- windows_all[i, ]
    res <- analyze_one_window(
      dat = dat,
      e_start = w$e_start,
      e_end = w$e_end,
      l_start = w$l_start,
      l_end = w$l_end,
      t_max = t_max
    )

    rt <- res$risk_tab
    ESS_E <- rt$ESS[rt$arm == "early"]
    ESS_L <- rt$ESS[rt$arm == "late"]

    data.frame(
      sim_setting = sim_name,
      family = w$family,
      window_id = w$window_id,
      window_name = w$window_name,
      e_start = w$e_start,
      e_end = w$e_end,
      l_start = w$l_start,
      l_end = w$l_end,
      OVL_L0 = res$OVL_L0,
      RR = res$RR,
      ESS_early = ESS_E,
      ESS_late = ESS_L
    )
  })

  all_rr_list[[sim_name]] <- bind_rows(res_each_sim)
}

RR_all <- bind_rows(all_rr_list) |>
  separate(sim_setting, into = c("Sim", "Scenario"), sep = "_") |>
  mutate(
    family = factor(
      family,
      levels = c("Sliding", "Late_disjoint", "Late_nested", "Gap_expansion", "Overlap_crossing")
    )
  ) |>
  group_by(family) |>
  arrange(window_id, .by_group = TRUE) |>
  mutate(
    window_label = paste0("[", e_start, "-", e_end, "] vs [", l_start, "-", l_end, "]"),
    window_label = factor(window_label, levels = unique(window_label))
  ) |>
  ungroup() |>
  mutate(
    SimScenario = factor(
      paste0(Sim, "_", Scenario),
      levels = c(
        "Sim0_A", "Sim0_B",
        "Sim1_A", "Sim1_B",
        "Sim2_A", "Sim2_B",
        "Sim3_A", "Sim3_B"
      )
    ),
    ESS_min = pmin(ESS_early, ESS_late)
  )

# Figure 3 panels

shape_vals <- c(
  "Sim0_A" = 1,
  "Sim0_B" = 2,
  "Sim1_A" = 16,
  "Sim1_B" = 17,
  "Sim2_A" = 15,
  "Sim2_B" = 18,
  "Sim3_A" = 0,
  "Sim3_B" = 5
)

linetype_vals <- c(
  "Sim0_A" = "solid",
  "Sim0_B" = "dashed",
  "Sim1_A" = "solid",
  "Sim1_B" = "dashed",
  "Sim2_A" = "dotdash",
  "Sim2_B" = "twodash",
  "Sim3_A" = "longdash",
  "Sim3_B" = "dotted"
)

pd <- position_dodge(width = 0.4)

p_rr <- ggplot(
  RR_all,
  aes(
    x = window_label,
    y = RR,
    group = SimScenario,
    color = SimScenario,
    shape = SimScenario,
    linetype = SimScenario
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(position = pd) +
  geom_point(position = pd, size = 2) +
  facet_wrap(~family, scales = "free_x") +
  labs(
    x = "Window setting",
    y = "12-month Risk Ratio (Early / Late)",
    color = "Simulation\nsetting",
    shape = "Simulation\nsetting",
    linetype = "Simulation\nsetting",
    title = "Primary analysis: Early-Late risk ratios across window families"
  ) +
  scale_shape_manual(values = shape_vals) +
  scale_linetype_manual(values = linetype_vals) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ovl <- ggplot(
  RR_all,
  aes(
    x = window_label,
    y = OVL_L0,
    group = SimScenario,
    color = SimScenario,
    shape = SimScenario,
    linetype = SimScenario
  )
) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 2) +
  facet_wrap(~family, scales = "free_x") +
  labs(
    x = "Window setting",
    y = "OVL(L0) between Early and Late",
    color = "Simulation\nsetting",
    shape = "Simulation\nsetting",
    linetype = "Simulation\nsetting",
    title = "Primary analysis: baseline covariate-support overlap across window families"
  ) +
  scale_shape_manual(values = shape_vals) +
  scale_linetype_manual(values = linetype_vals) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ess <- ggplot(
  RR_all,
  aes(
    x = window_label,
    y = ESS_min,
    group = SimScenario,
    color = SimScenario,
    shape = SimScenario,
    linetype = SimScenario
  )
) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 2) +
  facet_wrap(~family, scales = "free_x") +
  labs(
    x = "Window setting",
    y = "Minimum ESS (Early / Late)",
    color = "Simulation\nsetting",
    shape = "Simulation\nsetting",
    linetype = "Simulation\nsetting",
    title = "Primary analysis: effective sample size across window families"
  ) +
  scale_shape_manual(values = shape_vals) +
  scale_linetype_manual(values = linetype_vals) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_rr
p_ovl
p_ess

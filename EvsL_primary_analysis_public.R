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

############################################################
## Sim3 Monte Carlo bias analysis: Sliding windows only
############################################################

required_functions <- c(
  "simulate_dataset",
  "analyze_one_window",
  "make_windows_sliding"
)

missing_functions <- required_functions[!vapply(required_functions, exists, logical(1), mode = "function")]
if (length(missing_functions) > 0) {
  stop(
    "Missing required function(s): ",
    paste(missing_functions, collapse = ", "),
    "."
  )
}

required_objects <- c("t_max", "S_min", "S_max")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects) > 0) {
  stop(
    "Missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    "."
  )
}

windows_sliding_mc <- make_windows_sliding() %>%
  dplyr::mutate(
    family = "Sliding",
    window_name = paste0("Sliding_", window_id)
  ) %>%
  dplyr::arrange(window_id)

if (nrow(windows_sliding_mc) == 0) {
  stop("No sliding windows were generated by make_windows_sliding().")
}

run_one_rep_sim3_sliding <- function(n_sample,
                                     t_max,
                                     S_min,
                                     S_max) {
  dat_sim3_a <- simulate_dataset(
    n = n_sample,
    sim_type = 3,
    scenario = "A",
    t_max = t_max,
    S_min = S_min,
    S_max = S_max
  )

  dat_sim3_b <- simulate_dataset(
    n = n_sample,
    sim_type = 3,
    scenario = "B",
    t_max = t_max,
    S_min = S_min,
    S_max = S_max
  )

  scenario_data <- list(
    A = dat_sim3_a,
    B = dat_sim3_b
  )

  res_list <- vector("list", length = 2L * nrow(windows_sliding_mc))
  idx <- 1L

  for (scenario in c("A", "B")) {
    dat <- scenario_data[[scenario]]

    for (i in seq_len(nrow(windows_sliding_mc))) {
      w <- windows_sliding_mc[i, ]

      out <- tryCatch(
        analyze_one_window(
          dat = dat,
          e_start = w$e_start,
          e_end = w$e_end,
          l_start = w$l_start,
          l_end = w$l_end,
          t_max = t_max
        ),
        error = function(e) {
          warning(
            paste0(
              "analyze_one_window() failed for Sim3, Scenario ",
              scenario,
              ", window_id = ",
              w$window_id,
              ". Returning NA. Message: ",
              conditionMessage(e)
            ),
            call. = FALSE
          )
          list(RR = NA_real_)
        }
      )

      rr_hat <- suppressWarnings(as.numeric(out$RR)[1])

      if (length(rr_hat) == 0L || !is.finite(rr_hat)) {
        rr_hat <- NA_real_
      }

      res_list[[idx]] <- data.frame(
        Sim = "Sim3",
        Scenario = scenario,
        family = "Sliding",
        window_id = w$window_id,
        e_start = w$e_start,
        e_end = w$e_end,
        l_start = w$l_start,
        l_end = w$l_end,
        RR_hat = rr_hat,
        stringsAsFactors = FALSE
      )

      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(res_list)
}

mc_sim3_sliding <- function(n_rep = 1000,
                            n_sample = 5000,
                            t_max = t_max,
                            S_min = S_min,
                            S_max = S_max,
                            seed = 12345,
                            progress_every = 25) {
  if (!is.numeric(n_rep) || length(n_rep) != 1L || is.na(n_rep) || n_rep < 1) {
    stop("n_rep must be a single positive number.")
  }
  if (!is.numeric(n_sample) || length(n_sample) != 1L || is.na(n_sample) || n_sample < 1) {
    stop("n_sample must be a single positive number.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("seed must be a single numeric value.")
  }
  if (!is.numeric(progress_every) || length(progress_every) != 1L || is.na(progress_every) || progress_every < 1) {
    stop("progress_every must be a single positive number.")
  }

  n_rep <- as.integer(n_rep)
  n_sample <- as.integer(n_sample)
  progress_every <- as.integer(progress_every)

  set.seed(seed)

  all_res <- vector("list", length = n_rep)

  for (r in seq_len(n_rep)) {
    if (r == 1L || r %% progress_every == 0L || r == n_rep) {
      message("Sim3 Monte Carlo replication ", r, " / ", n_rep)
    }

    one_res <- run_one_rep_sim3_sliding(
      n_sample = n_sample,
      t_max = t_max,
      S_min = S_min,
      S_max = S_max
    ) %>%
      dplyr::mutate(
        rep = r,
        bias = RR_hat - 1
      )

    all_res[[r]] <- one_res
  }

  res_long <- dplyr::bind_rows(all_res)

  mc_summary <- res_long %>%
    dplyr::group_by(
      Sim, Scenario, family, window_id,
      e_start, e_end, l_start, l_end
    ) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_nonmissing = sum(!is.na(RR_hat)),
      RR_mean = mean(RR_hat, na.rm = TRUE),
      RR_sd = sd(RR_hat, na.rm = TRUE),
      bias_mean = mean(bias, na.rm = TRUE),
      bias_sd = sd(bias, na.rm = TRUE),
      RMSE = sqrt(mean(bias^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      RR_mean = ifelse(is.nan(RR_mean), NA_real_, RR_mean),
      RR_sd = ifelse(is.nan(RR_sd), NA_real_, RR_sd),
      bias_mean = ifelse(is.nan(bias_mean), NA_real_, bias_mean),
      bias_sd = ifelse(is.nan(bias_sd), NA_real_, bias_sd),
      RMSE = ifelse(is.nan(RMSE), NA_real_, RMSE),
      window_label = paste0(
        "[", e_start, "-", e_end, "] vs [", l_start, "-", l_end, "]"
      ),
      x_id = window_id,
      Scenario = factor(Scenario, levels = c("A", "B"))
    ) %>%
    dplyr::arrange(Scenario, window_id)

  list(
    draws = res_long,
    summary = mc_summary
  )
}

mc_res_sim3 <- mc_sim3_sliding(
  n_rep = 1000,
  n_sample = 5000,
  t_max = t_max,
  S_min = S_min,
  S_max = S_max,
  seed = 12345,
  progress_every = 25
)

mc_res_sim3$summary

x_breaks_sim3 <- sort(unique(mc_res_sim3$summary$x_id))
x_labels_sim3 <- mc_res_sim3$summary %>%
  dplyr::filter(Scenario == "A") %>%
  dplyr::arrange(window_id) %>%
  dplyr::pull(window_label)

base_text_size <- 14
base_title_size <- 16
facet_title_size <- 15

theme_base_mc <- theme_bw(base_size = base_text_size) +
  theme(
    plot.title = element_text(size = base_title_size, face = "bold"),
    axis.title.x = element_text(size = base_text_size + 1),
    axis.title.y = element_text(size = base_text_size + 1),
    axis.text.x = element_text(size = base_text_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = base_text_size),
    legend.title = element_text(size = base_text_size + 1, face = "bold"),
    legend.text = element_text(size = base_text_size),
    strip.text = element_text(size = facet_title_size, face = "bold"),
    panel.grid.major = element_line(linewidth = 0.4),
    panel.grid.minor = element_line(linewidth = 0.2)
  )

p_bias_sim3_ribbon <- ggplot(
  mc_res_sim3$summary,
  aes(
    x = x_id,
    y = bias_mean,
    group = Scenario,
    color = Scenario,
    linetype = Scenario
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(
    aes(
      ymin = bias_mean - bias_sd,
      ymax = bias_mean + bias_sd,
      fill = Scenario
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(shape = Scenario), size = 3, na.rm = TRUE) +
  scale_x_continuous(
    breaks = x_breaks_sim3,
    labels = x_labels_sim3
  ) +
  scale_color_manual(
    values = c(
      "A" = "#984EA3",
      "B" = "#F781BF"
    )
  ) +
  scale_fill_manual(
    values = c(
      "A" = "#984EA3",
      "B" = "#F781BF"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "A" = "longdash",
      "B" = "dotted"
    )
  ) +
  scale_shape_manual(
    values = c(
      "A" = 7,
      "B" = 8
    )
  ) +
  labs(
    x = "Sliding window (Early vs Late)",
    y = "MC mean bias of RR (RR_hat - 1)",
    color = "Scenario",
    fill = "Scenario",
    linetype = "Scenario",
    shape = "Scenario",
    title = "Sim3 (true RR = 1): MC bias \u00b1 SD across Sliding windows"
  ) +
  theme_base_mc

p_rmse_sim3 <- ggplot(
  mc_res_sim3$summary,
  aes(
    x = x_id,
    y = RMSE,
    group = Scenario,
    color = Scenario,
    shape = Scenario,
    linetype = Scenario
  )
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 3, na.rm = TRUE) +
  scale_x_continuous(
    breaks = x_breaks_sim3,
    labels = x_labels_sim3
  ) +
  scale_color_manual(
    values = c(
      "A" = "#984EA3",
      "B" = "#F781BF"
    )
  ) +
  scale_shape_manual(
    values = c(
      "A" = 7,
      "B" = 8
    )
  ) +
  scale_linetype_manual(
    values = c(
      "A" = "longdash",
      "B" = "dotted"
    )
  ) +
  labs(
    x = "Sliding window (Early vs Late)",
    y = "RMSE of RR (relative to 1)",
    color = "Scenario",
    shape = "Scenario",
    linetype = "Scenario",
    title = "Sim3 (true RR = 1): RMSE across Sliding windows"
  ) +
  theme_base_mc

p_bias_sim3_ribbon
p_rmse_sim3

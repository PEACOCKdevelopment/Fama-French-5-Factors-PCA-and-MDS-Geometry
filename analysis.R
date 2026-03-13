knitr::opts_chunk$set(
  echo = TRUE,        # <-- SHOW CODE (RPubs style)
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,    # <-- code + output together
  comment = "#>",     # <-- output prefix like in many RPubs examples
  fig.width = 8,
  fig.height = 5,
  fig.align = "center"
)
set.seed(123)

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_pkgs <- c(
  "tidyverse", "lubridate", "gridExtra", "smacof", "vegan", "psych", "tibble"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

ff5_raw <- readr::read_csv("FF5_monthly.csv", show_col_types = FALSE)

names(ff5_raw)
head(ff5_raw, 5)

ff5 <- ff5_raw %>%
  rename(DATE = 1) %>%
  mutate(DATE = as.character(DATE)) %>%
  filter(str_detect(DATE, "^[0-9]{6}$")) %>%
  mutate(date = ymd(paste0(DATE, "01"))) %>%
  select(date, `Mkt-RF`, SMB, HML, RMW, CMA) %>%
  arrange(date) %>%
  mutate(across(-date, as.numeric))

glimpse(ff5)
range(ff5$date)
summary(ff5)

diag_tbl <- tibble(
  start_date = min(ff5$date),
  end_date   = max(ff5$date),
  n_months   = nrow(ff5),
  n_factors  = ncol(ff5) - 1
)
knitr::kable(diag_tbl, caption = "Dataset diagnostics")

X <- ff5 %>% select(`Mkt-RF`, SMB, HML, RMW, CMA)

corr_mat <- cor(X, use = "pairwise.complete.obs")
round(corr_mat, 2)

knitr::kable(round(corr_mat, 2), caption = "Correlation matrix")

corr_df <- as.data.frame(corr_mat) %>%
  rownames_to_column("factor1") %>%
  pivot_longer(-factor1, names_to = "factor2", values_to = "corr")

ggplot(corr_df, aes(factor1, factor2, fill = corr)) +
  geom_tile() +
  geom_text(aes(label = round(corr, 2)), size = 4) +
  scale_fill_gradient2(limits = c(-1, 1)) +
  coord_equal() +
  labs(title = "Correlation heatmap", fill = "Corr") +
  theme_minimal()

corr_with_mkt <- sort(corr_mat[, "Mkt-RF"], decreasing = TRUE)
round(corr_with_mkt, 2)

X_scaled <- scale(X)
pca <- prcomp(X_scaled)

eig <- pca$sdev^2
prop <- eig / sum(eig)

pca_var_tbl <- tibble(
  PC = paste0("PC", 1:5),
  eigenvalue = eig,
  prop_var = prop,
  cum_var = cumsum(prop)
)

pca_var_tbl

knitr::kable(
  pca_var_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
  caption = "PCA: eigenvalues and explained variance"
)

scree_tbl <- pca_var_tbl %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:5)))

ggplot(scree_tbl, aes(PC, prop_var, group = 1)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_line(color = "grey30") +
  geom_point(color = "grey30", size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Scree plot: variance explained by PCs",
    x = NULL,
    y = "Explained variance"
  ) +
  theme_minimal()

loadings <- pca$rotation[, 1:3] %>%
  as.data.frame() %>%
  rownames_to_column("Factor")

knitr::kable(
  loadings %>% mutate(across(-Factor, ~ round(.x, 3))),
  caption = "PCA loadings (PC1-PC3)"
)

contrib_tbl <- as_tibble(pca$rotation[, 1:3], rownames = "Factor") %>%
  pivot_longer(-Factor, names_to = "PC", values_to = "loading") %>%
  group_by(PC) %>%
  mutate(contrib = 100 * (loading^2) / sum(loading^2)) %>%
  ungroup()

plot_contrib <- function(pc_name) {
  ggplot(
    contrib_tbl %>% filter(PC == pc_name) %>% arrange(desc(contrib)),
    aes(x = reorder(Factor, contrib), y = contrib)
  ) +
    geom_col(fill = "steelblue", alpha = 0.85) +
    coord_flip() +
    labs(title = paste("Contribution to", pc_name), x = NULL, y = "Contribution (%)") +
    theme_minimal()
}

gridExtra::grid.arrange(
  plot_contrib("PC1"),
  plot_contrib("PC2"),
  plot_contrib("PC3"),
  ncol = 1
)

scores <- as_tibble(pca$x) %>%
  mutate(date = ff5$date)

ggplot(scores, aes(PC1, PC2)) +
  geom_point(alpha = 0.5) +
  labs(title = "Months in PC space (PC1 vs PC2)") +
  theme_minimal()

scores_12 <- as_tibble(pca$x[, 1:2], .name_repair = "minimal") %>%
  rename(PC1 = 1, PC2 = 2)

loadings_12 <- as_tibble(pca$rotation[, 1:2], rownames = "Factor")

arrow_scale <- min(
  diff(range(scores_12$PC1)) / diff(range(loadings_12$PC1)),
  diff(range(scores_12$PC2)) / diff(range(loadings_12$PC2))
) * 0.25

loadings_plot <- loadings_12 %>%
  mutate(
    PC1 = PC1 * arrow_scale,
    PC2 = PC2 * arrow_scale
  )

ggplot(scores_12, aes(PC1, PC2)) +
  geom_point(alpha = 0.5, color = "grey40") +
  geom_segment(
    data = loadings_plot,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    color = "steelblue",
    arrow = arrow(length = grid::unit(0.2, "cm"))
  ) +
  geom_text(
    data = loadings_plot,
    aes(x = PC1, y = PC2, label = Factor),
    inherit.aes = FALSE,
    color = "steelblue",
    nudge_y = 0.05,
    size = 4
  ) +
  labs(title = "PCA biplot: PC1 vs PC2") +
  theme_minimal()

dist_mat <- dist(scale(X))
mds_cmd <- cmdscale(dist_mat, k = 2)

mds_df <- tibble(
  Dim1 = mds_cmd[, 1],
  Dim2 = mds_cmd[, 2]
)

ggplot(mds_df, aes(Dim1, Dim2)) +
  geom_point(alpha = 0.5, color = "grey40") +
  labs(title = "Classical MDS (2D)") +
  theme_minimal()

mds_smacof <- smacof::mds(dist_mat, ndim = 2, type = "ratio")
mds_smacof$stress

knitr::kable(
  tibble(stress_2D = round(mds_smacof$stress, 4)),
  caption = "SMACOF MDS stress (2D)"
)

plot(mds_smacof, main = "SMACOF MDS (stress-based)")

stress_curve <- tibble(
  ndim = 1:5,
  stress = map_dbl(1:5, ~ smacof::mds(dist_mat, ndim = .x, type = "ratio")$stress)
)

knitr::kable(
  stress_curve %>% mutate(stress = round(stress, 4)),
  caption = "Stress vs number of dimensions"
)

ggplot(stress_curve, aes(ndim, stress)) +
  geom_line() +
  geom_point(size = 2) +
  labs(title = "SMACOF MDS: stress curve") +
  theme_minimal()

d_orig <- as.vector(dist(scale(X)))
d_pca  <- as.vector(dist(as.matrix(scores[, c("PC1","PC2")])))
d_mds  <- as.vector(dist(mds_cmd))

cor(d_orig, d_pca)
cor(d_orig, d_mds)

knitr::kable(
  tibble(
    comparison = c("Original vs PCA(2D)", "Original vs MDS(2D)"),
    correlation = round(c(cor(d_orig, d_pca), cor(d_orig, d_mds)), 4)
  ),
  caption = "Distance preservation"
)

proc <- vegan::procrustes(
  as.matrix(scores[, c("PC1","PC2")]),
  as.matrix(mds_cmd),
  symmetric = TRUE
)

proc_overlay <- bind_rows(
  tibble(Dim1 = proc$X[, 1], Dim2 = proc$X[, 2], Space = "PCA (2D)"),
  tibble(Dim1 = proc$Yrot[, 1], Dim2 = proc$Yrot[, 2], Space = "MDS aligned")
)

ggplot(proc_overlay, aes(Dim1, Dim2, color = Space)) +
  geom_point(alpha = 0.35, size = 1) +
  labs(title = "Procrustes alignment: PCA vs MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()

summary(proc)

rot_pca <- psych::principal(X_scaled, nfactors = 3, rotate = "varimax")
print(rot_pca$loadings, cutoff = 0.3)

portfolio_lines <- readr::read_lines("25_Portfolios_5x5.csv")
ew_monthly_marker <- which(str_detect(portfolio_lines, "Average Equal Weighted Returns -- Monthly"))[1]

if (is.na(ew_monthly_marker) || ew_monthly_marker <= 3) {
  stop("Cannot locate the monthly value-weighted block in 25_Portfolios_5x5.csv")
}

vw_last_data_line <- max(which(str_detect(portfolio_lines[1:(ew_monthly_marker - 1)], "^[0-9]{6},")))

vw_monthly_text <- paste(portfolio_lines[1:vw_last_data_line], collapse = "\n")

port25_vw <- readr::read_csv(
  I(vw_monthly_text),
  show_col_types = FALSE,
  na = c("-99.99", "-999")
) %>%
  rename(DATE = 1) %>%
  mutate(
    DATE = as.character(DATE),
    date = ymd(paste0(DATE, "01"))
  ) %>%
  filter(str_detect(DATE, "^[0-9]{6}$")) %>%
  select(-DATE)

rf_monthly <- ff5_raw %>%
  rename(DATE = 1) %>%
  mutate(DATE = as.character(DATE)) %>%
  filter(str_detect(DATE, "^[0-9]{6}$")) %>%
  mutate(date = ymd(paste0(DATE, "01"))) %>%
  select(date, RF) %>%
  mutate(RF = as.numeric(RF))

pc_scores <- scores %>% select(date, PC1, PC2, PC3)
ff5_original <- ff5 %>% select(date, `Mkt-RF`, SMB, HML, RMW, CMA)

analysis_panel <- port25_vw %>%
  inner_join(pc_scores, by = "date") %>%
  inner_join(ff5_original, by = "date") %>%
  inner_join(rf_monthly, by = "date")

portfolio_names <- setdiff(names(port25_vw), "date")

section9_diag <- tibble(
  start_date = min(analysis_panel$date),
  end_date = max(analysis_panel$date),
  n_months = nrow(analysis_panel),
  n_portfolios = length(portfolio_names),
  predictors_PCA = 3,
  predictors_FF5 = 5
)

knitr::kable(section9_diag, caption = "Chapter 9 data panel diagnostics")

comparison_results <- map_dfr(portfolio_names, function(pname) {
  model_data <- analysis_panel %>%
    transmute(
      ret_excess = .data[[pname]] - RF,
      PC1, PC2, PC3,
      `Mkt-RF`, SMB, HML, RMW, CMA
    ) %>%
    drop_na()

  fit_pca <- lm(ret_excess ~ PC1 + PC2 + PC3, data = model_data)
  fit_ff5 <- lm(ret_excess ~ `Mkt-RF` + SMB + HML + RMW + CMA, data = model_data)

  co_pca <- coef(fit_pca)
  co_ff5 <- coef(fit_ff5)

  s_pca <- summary(fit_pca)
  s_ff5 <- summary(fit_ff5)

  tibble(
    portfolio = pname,
    alpha_PCA = unname(co_pca[1]),
    beta_PC1 = unname(co_pca[2]),
    beta_PC2 = unname(co_pca[3]),
    beta_PC3 = unname(co_pca[4]),
    alpha_FF5 = unname(co_ff5[1]),
    beta_MKT = unname(co_ff5[2]),
    beta_SMB = unname(co_ff5[3]),
    beta_HML = unname(co_ff5[4]),
    beta_RMW = unname(co_ff5[5]),
    beta_CMA = unname(co_ff5[6]),
    R2_PCA = s_pca$r.squared,
    adjR2_PCA = s_pca$adj.r.squared,
    R2_FF5 = s_ff5$r.squared,
    adjR2_FF5 = s_ff5$adj.r.squared,
    delta_adjR2 = s_pca$adj.r.squared - s_ff5$adj.r.squared,
    R2_retention = if_else(s_ff5$r.squared > 0, s_pca$r.squared / s_ff5$r.squared, NA_real_),
    n_obs = nobs(fit_pca)
  )
})

knitr::kable(
  comparison_results %>%
    select(portfolio, beta_PC1, beta_PC2, beta_PC3, R2_PCA, adjR2_PCA, R2_FF5, adjR2_FF5, delta_adjR2) %>%
    mutate(across(-portfolio, ~ round(.x, 4))),
  caption = "Portfolio-level comparison: PCA latent model vs original FF5"
)

model_compare_summary <- comparison_results %>%
  summarize(
    mean_R2_PCA = mean(R2_PCA, na.rm = TRUE),
    mean_R2_FF5 = mean(R2_FF5, na.rm = TRUE),
    mean_adjR2_PCA = mean(adjR2_PCA, na.rm = TRUE),
    mean_adjR2_FF5 = mean(adjR2_FF5, na.rm = TRUE),
    mean_R2_retention = mean(R2_retention, na.rm = TRUE),
    median_R2_retention = median(R2_retention, na.rm = TRUE),
    mean_delta_adjR2 = mean(delta_adjR2, na.rm = TRUE),
    share_adjR2_within_95pct = mean(adjR2_PCA >= 0.95 * adjR2_FF5, na.rm = TRUE),
    share_adjR2_PCA_gt_FF5 = mean(adjR2_PCA > adjR2_FF5, na.rm = TRUE)
  )

knitr::kable(
  model_compare_summary %>% mutate(across(everything(), ~ round(.x, 4))),
  caption = "Chapter 9 effectiveness summary"
)

ggplot(comparison_results, aes(adjR2_FF5, adjR2_PCA)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2.2, alpha = 0.85, color = "steelblue") +
  labs(
    title = "Adjusted R-squared: PCA latent model vs original FF5",
    x = "Adjusted R-squared (FF5)",
    y = "Adjusted R-squared (PCA latent factors)"
  ) +
  theme_minimal()

beta_long <- comparison_results %>%
  select(portfolio, beta_PC1, beta_PC2, beta_PC3) %>%
  pivot_longer(-portfolio, names_to = "component", values_to = "beta") %>%
  mutate(component = str_remove(component, "beta_"))

portfolio_grid <- tibble(portfolio = unique(beta_long$portfolio)) %>%
  mutate(
    size = case_when(
      portfolio == "SMALL LoBM" ~ 1L,
      portfolio == "SMALL HiBM" ~ 1L,
      portfolio == "BIG LoBM" ~ 5L,
      portfolio == "BIG HiBM" ~ 5L,
      TRUE ~ as.integer(str_match(portfolio, "^ME([1-5]) BM([1-5])$")[, 2])
    ),
    bm = case_when(
      portfolio == "SMALL LoBM" ~ 1L,
      portfolio == "SMALL HiBM" ~ 5L,
      portfolio == "BIG LoBM" ~ 1L,
      portfolio == "BIG HiBM" ~ 5L,
      TRUE ~ as.integer(str_match(portfolio, "^ME([1-5]) BM([1-5])$")[, 3])
    ),
    size_label = factor(size, levels = 5:1, labels = c("Big", "4", "3", "2", "Small")),
    bm_label = factor(bm, levels = 1:5, labels = c("LoBM", "2", "3", "4", "HiBM"))
  )

beta_heatmap <- beta_long %>%
  left_join(portfolio_grid, by = "portfolio")

ggplot(beta_heatmap, aes(bm_label, size_label, fill = beta)) +
  geom_tile(color = "white") +
  facet_wrap(~ component, nrow = 1) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c") +
  labs(
    title = "PCA latent-factor exposures across the 5x5 size-value grid",
    x = "Book-to-Market bucket",
    y = "Size bucket",
    fill = "Beta"
  ) +
  theme_minimal()

pkg_versions <- tibble(
  package = c("R", required_pkgs),
  version = c(
    as.character(getRversion()),
    vapply(required_pkgs, function(p) as.character(packageVersion(p)), character(1))
  )
)

knitr::kable(pkg_versions, caption = "R and package versions used in this run")
sessionInfo()

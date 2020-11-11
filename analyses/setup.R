
# load packages -----------------------------------------------------------

library(tidyverse)
library(rstatix)
library(here)
source(here::here("analyses/WRS-v35.txt"))

aov_2levels <- function(df, factor, model, level1, level2) {
  anova_table <- df %>% 
    dplyr::filter(!!rlang::sym(factor) %in% c(level1, level2)) %>% 
    rstatix::anova_test(model) %>% 
    rstatix::get_anova_table(correction = "auto")
  
  df_w <- df %>% 
    tidyr::pivot_wider(names_from = !!rlang::sym(factor), values_from = attr(anova_table, "args")$dv)
  g <- MBESS::smd(df_w[[level1]], df_w[[level2]], Unbiased = T)
  
  res <- anova_table %>% 
    as_tibble() %>% 
    dplyr::filter(Effect == factor) %>% 
    dplyr::mutate(level1 = level1, level2 = level2, t = sqrt(`F`), g = abs(g))
  
  return(res)
}

multi_aov <- function(anova_table, factor) {
  lvls <- attr(anova_table, "args")$data[[factor]] %>% levels
  contrasts <- combn(lvls, m = 2)  %>% 
    split(col(.))
  
  df <- attr(anova_table, "args")$data
  model <- attr(anova_table, "args")$formula
  
  res <- contrasts %>% 
    purrr::map(~aov_2levels(df, factor, model, .x[1], .x[2]))%>% 
    dplyr::bind_rows() %>% 
    dplyr::select(level1, level2, DFn:p, ges, t, g) %>% 
    rstatix::adjust_pvalue(p.col = "p", method = "holm")

  return(res)
}

format_pval <- function(pval) {
  sapply(pval, function(x) {
    if (is.na(x)) {
      return(NA_character_)
    }
    if (x < .001) {
      pval_txt <- "_p_ < .001"
    }
    else if (x == 1) {
      pval_txt <- "_p_ = 1.00"
    }
    else {
      pval_new <- 
        sprintf("%.3f", x) %>% 
        stringr::str_sub(start = 2)
      pval_txt <- sprintf("_p_ = %s", pval_new)
    }
    return(pval_txt)
  })
}

format_aov <- function(aov_tbl, grouped = F, p_adj = F) {
  pval <- "{format_pval(p)}, "
  if (p_adj) pval <- "{format_pval(p.adj)}, "

  txt <- aov_tbl %>% 
    stringr::str_glue_data("_F_({DFn}, {DFd}) = {sprintf('%.2f', F)}, ", pval, "$\\eta^2_G$ = {sprintf('%.2f', ges)}")
  
  if (!grouped) {
    names(txt) <- aov_tbl$Effect
    return(as.list(txt))
  }
  names(txt) <- 
    str_glue("{Effect}@{grp}",
             Effect = aov_tbl$Effect,
             grp = dplyr::pull(aov_tbl, 1))
  return(as.list(txt))
}

format_mc <- function(mc_tbl, es = F) {
  pattern <- "_t_({DFd}) = {sprintf('%.2f', abs(t))}, {format_pval(p.adj)}"
  if (es) {
    pattern <- "_t_({DFd}) = {sprintf('%.2f', abs(t))}, {format_pval(p.adj)}, g = {sprintf('%.2f', abs(g))}"
  }

  txt <- mc_tbl %>% 
    stringr::str_glue_data(pattern)

  names(txt) <- str_glue_data(mc_tbl, "{level1}^{level2}")
  
  return(as.list(txt))
}

format_t <- function(t_tbl, p_adj = T, grouped = F, es = F) {
  pval <- "{format_pval(p)}"
  if (p_adj) {
    pval <- "{format_pval(p.adj)}"
  }
  
  format_df <- function(x) {
    if (is.double(x)) {
      return(sprintf("%.2f", x))
    }
    return(x)
  }

  txt <- t_tbl %>% 
    stringr::str_glue_data("_t_({format_df(df)}) = {sprintf('%.2f', abs(statistic))}, ", pval)
  if (es) {
    txt <- t_tbl %>% 
      stringr::str_glue_data("_t_({format_df(df)}) = {sprintf('%.2f', abs(statistic))}, ", 
                             pval, 
                             ", g = {sprintf('%.2f', abs(g))}")
  }
  
  if (!grouped) {
    names(txt) <- str_glue_data(t_tbl, "{group1}^{group2}")
    return(as.list(txt))
  }
  names(txt) <- 
    str_glue("{level1}^{level2}@{grp}",
             level1 = t_tbl$group1,
             level2 = t_tbl$group2,
             grp = dplyr::pull(t_tbl, 1))
  
  return(as.list(txt))
}

# for creating summary tables
mean_cl <- function(x){
  res <- Hmisc::smean.cl.normal(x)
  mean <- sprintf("%.2f", res["Mean"])
  ci <- sprintf("[%.2f, %.2f]", res['Lower'], res['Upper'])
  return(c(mean, ci))
}

# load data sets ----------------------------------------------------------

df_e1 <- fs::dir_ls(path = here::here("exp/data/exp1"), regexp = "[[:digit:]]+\\.csv") %>% 
  purrr::map_dfr(readr::read_csv) %>% 
  dplyr::mutate(suppression = if_else(participant %% 2 == 0, "TS", "DS"))

df_ques_e1 <- fs::dir_ls(path = here::here("exp/data/exp1"), regexp = "[[:digit:]]+_ques\\.csv") %>% 
  purrr::map_dfr(readr::read_csv) %>% 
  dplyr::mutate(suppression = if_else(participant %% 2 == 0, "TS", "DS"))

df_e2 <- fs::dir_ls(path = here::here("exp/data/exp2"), regexp = "[[:digit:]]+\\.csv") %>% 
  purrr::map_dfr(readr::read_csv) %>% 
  dplyr::filter(participant != 3) %>% 
  dplyr::mutate(suppression = if_else(participant %% 2 == 0, "TS", "DS"))

df_ques_e2 <- fs::dir_ls(path = here::here("exp/data/exp2"), regexp = "[[:digit:]]+_ques\\.csv") %>% 
  purrr::map_dfr(readr::read_csv) %>% 
  dplyr::filter(participant != 3) %>% 
  dplyr::mutate(suppression = if_else(participant %% 2 == 0, "TS", "DS")) %>% 
  dplyr::rename("type" = condition, "status" = surppress)

df_rcll_l_sub <- fs::dir_ls(path = here::here("exp/data/exp2"), regexp = "[[:digit:]]+_sRT\\.csv") %>% 
  str_subset(pattern = "/(3|12|27)_sRT\\.csv", negate = T) %>% 
  purrr::map_dfr(readr::read_csv) %>% 
  dplyr::mutate(suppression = if_else(participant %% 2 == 0, "TS", "DS"))

df_dp_raw_e2 <- purrr::map_dfr(
  .x = fs::dir_ls(here::here("exp/data/exp2"), glob = "*.xlsx"),
  .f = function(x) {
    par_id <- stringr::str_extract(x, "\\d+_TNT") %>% readr::parse_number()
    snms <- readxl::excel_sheets(x)
    rng <- "A2:J137"
    colnms <- c("ID", "below", "letter", "above", "congruency", "position", "order", "key", "corr", "RT")
    if (length(snms) > 10) {
      rng <- "A1:P137"
      colnms <- TRUE
    }
    dp_snm <- str_subset(snms, "(DP)|(dot)")
    df <- readxl::read_xlsx(path = x, sheet = dp_snm, range = rng, col_names = colnms) %>% 
      dplyr::mutate(participant = par_id)
    if (length(snms) > 10) {
      df <- df %>% 
        dplyr::select(participant, ID, below, letter, above, congruency, position, order,
                      "key" = key_resp_43.keys_raw, 
                      "corr" = key_resp_43.corr_mean,
                      "RT" = key_resp_43.rt_mean)
    }
    return(df)
  }
) %>% dplyr::left_join(
    y = df_e2 %>% dplyr::select(participant, item_id, status),
    by = c("participant", "ID" = "item_id")
  ) %>% 
  dplyr::mutate(status = replace_na(status, "filler"))


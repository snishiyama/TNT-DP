library(tidyverse)
library(rstatix)
library(here)
source(here::here("analyses/WRS-v35.txt"))

multi_aov <- function(df, factor, model, levels) {
  aov_2levels <- function(df, factor, model, level1, level2) {
    df %>% 
      dplyr::filter(!!rlang::sym(factor) %in% c(level1, level2)) %>% 
      rstatix::anova_test(model) %>% 
      rstatix::get_anova_table(correction = "auto") %>% 
      as_tibble() %>% 
      dplyr::filter(Effect == factor) %>% 
      dplyr::mutate(level1 = level1, level2 = level2)
  }

  contrasts <- combn(levels, m = 2)  %>% 
    split(col(.))

  res <- contrasts %>% 
    purrr::map(~aov_2levels(df, factor, model, .x[1], .x[2]))%>% 
    dplyr::bind_rows() %>% 
    dplyr::select(level1, level2, DFn:p, ges) %>% 
    rstatix::adjust_pvalue(p.col = "p", method = "holm") %>% 
    dplyr::mutate(t = sqrt(`F`))

  return(res)
}


df_e1 <- readr::read_csv(here::here("exp/data/sotsuron_exp1_part2.csv")) %>% 
  dplyr::select(-1)

df_e2 <- readr::read_csv(here::here("exp/data/TNT_exp2_sotsuron.csv")) %>% 
  dplyr::mutate(suppression = suppression - 1L) # Exp1の条件分けが0 / 1だから

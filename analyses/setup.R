library(tidyverse)
library(rstatix)
library(here)

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

df_recall_e1 <- df_e1 %>% 
  dplyr::select(participant, suppression, 
                think = recall_think, 
                nothink = recall_nothink, 
                baseline = recall_base, 
                substitute = recall_sub) %>% 
  tidyr::pivot_longer(think:substitute, names_to = "status", values_to = "rate")

# anova t vs nt vs b
model_recall_e1 <- rate ~ suppression*status + Error(participant/status)
aov_recall_e1 <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::anova_test(model_recall_e1) %>% 
  rstatix::get_anova_table(correction = "auto")


# multiple comparisons
multi_aov(df_recall_e1, "status", model_recall_e1, c("think", "nothink", "baseline"))

# multiple comparisons with emmeans
df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  dplyr::mutate(participant = as_factor(participant),
                suppression = as_factor(suppression),
                status = factor(status, levels = c("think", "nothink", "baseline"))) %>% 
  aov(rate ~ suppression*status+Error(participant/status), 
      contrasts = list(suppression = "contr.sum", status = "contr.sum"), # type 3 SS
      data = .) %>% 
  emmeans::emmeans(~ status) %>% 
  emmeans::contrast(method = "pairwise", adjust = "holm")
  

pwt_recall_e1 <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::pairwise_t_test(rate ~ status, paired = T, p.adjust.method = "holm", detailed = T)

rate <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  dplyr::pull(rate)
  
# source(here::here("analyses/anovakun_485.txt"))
# aovkun <- df_recall_e1 %>%
#   dplyr::filter(status != "substitute") %>%
#   anovakun("AsB", long = T, geta = T)

# calculate effect size
# pipeR::pipeline({
#   exp1_res_anova_rate$`POST ANALYSES`$B$bontab
#   ? .
#   split(.$pair)
#   purrr::map_if(.p = ~ .$adj.p < 0.05, .f = ~MBESS::ci.sm(ncp = .$t[[1]], N = .$df[[1]] + 1))
# })
# t.test on substitute
df_t_sub_recall_e1 <- df_recall_e1 %>% 
  dplyr::filter(status == "substitute") %>% 
  rstatix::t_test(rate ~ suppression)

# d <- ttest_sub_exp1$statistic[[1]] * sqrt((18+18)/(18*18))
# MBESS::ci.smd(smd = d, n.1 = 18, n.2 = 18)
MBESS::ci.smd(ncp = ttest_sub_exp1$statistic[[1]], n.1 = 18, n.2 = 18)

df_ques_e1 <- df_e1 %>% 
  dplyr::select(participant, suppression, starts_with("think"), starts_with("nothink")) %>% 
  tidyr::pivot_longer(cols = -c(suppression, participant), 
                      names_to = c("condition", "type"), 
                      names_sep = "_",
                      values_to = "rating")

model_ques_e1 <- rating ~ suppression*type+Error(participant/type)
l_aov_ques_e1 <- df_ques_e1 %>% 
  split(.$condition) %>% 
  purrr::map(.f = rstatix::anova_test, model_ques_e1) %>% 
  purrr::map(.f = rstatix::get_anova_table, correction = "auto")

mc_ques_think_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "think") %>% 
  multi_aov("type", model_ques_e1, c("cue", "tar", "sub"))

sme_ques_nt_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(rating ~ type + Error(participant/type)) %>% 
  rstatix::get_anova_table(correction = "auto")

pwt_ques_nt_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)

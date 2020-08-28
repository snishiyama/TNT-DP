library(tidyverse)
library(rstatix)
library(here)

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
aov_recall_e1 <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::anova_test(rate ~ suppression*status + Error(participant/status))  %>% 
  rstatix::get_anova_table(correction = "auto")

pwt_recall_e1 <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::pairwise_t_test(rate ~ status, paired = T, p.adjust.method = "holm", detailed = T)

rate <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  dplyr::pull(rate)

df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::emmeans_test(rate ~ status, p.adjust.method = "holm", detailed = T)
  
source(here::here("analyses/anovakun_485.txt"))
aovkun <- df_recall_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  anovakun("AsB", long = T, geta = T)

# 0.0898/sqrt(0.0050*(2/36))
# 
# 0.898 / 4.4901

# temp2 <- df_recall_e1 %>% 
#   dplyr::filter(!status %in% c("substitute", "think")) %>% 
#   tidyr::pivot_wider(names_from = status, values_from = rate) %>% 
#   dplyr::mutate(diff = nothink - baseline) %>% 
#   dplyr::summarise(mean = mean(diff), sd = sd(diff), se = sd/sqrt(36), t = mean/se)
# 
# temp <- df_recall_e1 %>% 
#   dplyr::filter(!status %in% c("substitute", "think")) %>% 
#   rstatix::anova_test(rate ~ status + Error(participant/status), detailed = T)
# 
# temp <- df_recall_e1 %>% 
#   dplyr::filter(!status %in% c("substitute", "think")) %>% 
#   dplyr::mutate(participant = as_factor(participant)) %>% 
#   lm(rate ~ status + participant, data = .) %>% 
#   car::Anova(type = 3)
# 
# bon.SE <- sqrt((bon.denomi[comb.frame[1, ]] + bon.denomi[comb.frame[2, ]]) * bon.Ve)
# bon.t <- abs(bon.delta / bon.SE)# ｔ値は絶対値を取る
# bon.p <- pt(bon.t, bon.df, lower.tail = FALSE) * 2# 両側確率
# cont.N <- table(dat[, 2:nchar(design)])# 各セルのデータ数
# bon.denomi <- apply(1/cont.N, bon.num, mean) / (prod(factlevels)/2)# セルごとに重み付けしたデータ数の平均を残りの条件数で割ったもの
# othink-baseline  -0.0760   4.3749  34  0.0001  0.0002  nothink < baseline *
  # (aovkun$`POST ANALYSES`$B$bontab$difference[2] / aovkun$`POST ANALYSES`$B$bontab$t[2])^2 / (temp$`Sum Sq`[4] / 35)
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

l_aov_ques_e1 <- df_ques_e1 %>% 
  split(.$condition) %>% 
  purrr::map(.f = rstatix::anova_test, rating ~ suppression*type+Error(participant/type)) %>% 
  purrr::map(.f = rstatix::get_anova_table, correction = "auto")

pwt_ques_think_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "think") %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)

sme_ques_nt_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(rating ~ type + Error(participant/type)) %>% 
  rstatix::get_anova_table(correction = "auto")

pwt_ques_nt_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)

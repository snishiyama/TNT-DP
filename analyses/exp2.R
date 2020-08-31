
# setup -------------------------------------------------------------------

if (!exists("df_e2")) {
  source(here::here("analyses/setup.R"))
}


# questionnaire -----------------------------------------------------------

df_ques_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, contains("rating")) %>% 
  tidyr::pivot_longer(cols = contains("rating"), 
                      names_to = c("condition", "type"), 
                      names_sep = "_rating.", 
                      values_to = "rating")

# ANOVA
model_ques_e2 <- rating ~ suppression*type+Error(participant/type)
l_aov_ques_e2 <- df_ques_e2 %>% 
  split(.$condition) %>% 
  purrr::map(.f = rstatix::anova_test, model_ques_e2) %>% 
  purrr::map(.f = rstatix::get_anova_table, correction = "auto")

# multiple comparisons for cue types of think items
mc_ques_think_e2 <- df_ques_e2 %>% 
  dplyr::filter(condition == "t") %>% 
  multi_aov("type", model_ques_e2, c("cue", "tar", "sub"))

df_es_think_ques_e2 <- df_ques_e2 %>% 
  dplyr::filter(condition == "t") %>% 
  tidyr::pivot_wider(names_from = type, values_from = rating) %>% 
  dplyr::summarise(cue_sub = MBESS::smd(cue, sub, Unbiased = T),
                   cue_tar = MBESS::smd(cue, tar, Unbiased = T),
                   tar_sub = MBESS::smd(tar, sub, Unbiased = T)) %>% 
  tidyr::pivot_longer(cols = cue_sub:tar_sub, names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# group-by simple main effect analyses for no-think items
sme_ques_nt_e2 <- df_ques_e2 %>% 
  dplyr::filter(condition == "nt") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(rating ~ type + Error(participant/type)) %>% 
  rstatix::get_anova_table(correction = "auto")

# pairwise t test 
pwt_ques_nt_e2 <- df_ques_e2 %>% 
  dplyr::filter(condition == "nt") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)

# effect sizes
df_es_nt_ques_e2 <- df_ques_e2 %>% 
  dplyr::filter(condition == "nt") %>% 
  tidyr::pivot_wider(names_from = type, values_from = rating) %>% 
  dplyr::group_by(suppression) %>% 
  dplyr::summarise(cue_sub = MBESS::smd(cue, sub, Unbiased = T),
                   cue_tar = MBESS::smd(cue, tar, Unbiased = T),
                   sub_tar = MBESS::smd(tar, sub, Unbiased = T)) %>% 
  tidyr::pivot_longer(cols = -suppression, names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# reporting
purrr::map(l_aov_ques_e2, format_aov)
dplyr::left_join(mc_ques_think_e2, df_es_think_ques_e2, by = c("level1", "level2")) %>% 
  format_mc(es = T)
dplyr::left_join(pwt_ques_nt_e2, df_es_nt_ques_e2, by = c("suppression", "group1" = "level1", "group2" = "level2")) %>% 
  format_pwt(grouped = T, es = T)


# recall rate -------------------------------------------------------------

df_rcll_r_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, ends_with("recall")) %>% 
  dplyr::mutate_at(vars(ends_with("recall")), ~ . / 100) %>% 
  dplyr::rename(think = t_recall, nothink = nt_recall,
                baseline = b_recall, substitute = s_recall) %>% 
  tidyr::pivot_longer(think:substitute, names_to = "status", values_to = "rate")

# anova t vs nt vs b
model_rcll_r_e2 <- rate ~ suppression*status + Error(participant/status)
aov_rcll_r_e2 <- df_rcll_r_e2 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::anova_test(model_rcll_r_e2) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_r_e2 <- multi_aov(df_rcll_r_e2, "status", model_rcll_r_e2, c("think", "nothink", "baseline"))

# Effect size
df_es_rcll_r_e2 <- df_rcll_r_e2 %>% 
  tidyr::pivot_wider(names_from = status, values_from = rate) %>% 
  dplyr::summarise(think_nothink = MBESS::smd(think, nothink, Unbiased = T),
                   think_baseline = MBESS::smd(think, baseline, Unbiased = T),
                   nothink_baseline = MBESS::smd(nothink, baseline, Unbiased = T)) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# t.test on substitute
ttest_sub_rcll_r_e2 <- df_rcll_r_e2 %>% 
  dplyr::filter(status == "substitute") %>% 
  rstatix::t_test(rate ~ suppression, paired = F, var.equal = F)

l_sub_rcll_r_e2 <- df_rcll_r_e2 %>% 
  dplyr::filter(status == "substitute") %>% 
  split(.$suppression)
g_rcll_r_sub_e2 <- MBESS::smd(l_sub_rcll_r_e2$`0`$rate, l_sub_rcll_r_e2$`1`$rate)


# reporting
format_aov(aov_rcll_r_e2)
dplyr::left_join(mc_rcll_r_e2, df_es_rcll_r_e2, by = c("level1", "level2")) %>% 
  format_mc(es = T)
ttest_sub_rcll_r_e2["g"] <- g_rcll_r_sub_e2#, df_es_nt_ques_e1, by = c("suppression", "group1" = "level1", "group2" = "level2")) %>% 
format_t(ttest_sub_rcll_r_e2, es = T)

# recall latency ----------------------------------------------------------

df_rcll_d_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, ends_with("delay")) %>% 
  dplyr::rename(think = t_delay, nothink = nt_delay, baseline = b_delay) %>% 
  tidyr::pivot_longer(think:baseline, names_to = "status", values_to = "delay") %>% 
  tidyr::drop_na() # remove data of id 27

# anova t vs nt vs b
model_rcll_d_e2 <- delay ~ suppression*status + Error(participant/status)
aov_rcll_d_e2 <- df_rcll_d_e2 %>% 
  rstatix::anova_test(model_rcll_d_e2) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_d_e2 <- multi_aov(df_rcll_d_e2, "status", model_rcll_d_e2, c("think", "nothink", "baseline"))

# Effect size
df_es_rcll_d_e2 <- df_rcll_d_e2 %>% 
  tidyr::pivot_wider(names_from = status, values_from = delay) %>% 
  dplyr::summarise(think_nothink = MBESS::smd(think, nothink, Unbiased = T),
                   think_baseline = MBESS::smd(think, baseline, Unbiased = T),
                   nothink_baseline = MBESS::smd(nothink, baseline, Unbiased = T)) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# group-by simple main effect analyses for no-think items
sme_rcll_d_e2 <- df_rcll_d_e2 %>% 
  dplyr::group_by(status) %>% 
  rstatix::anova_test(delay ~ suppression) %>% 
  rstatix::get_anova_table(correction = "auto")

# substitute
df_rcll_l_sub <- df_e2 %>% 
  dplyr::select(participant, suppression, s_RT)
ttest_sub_rcll_l_e2 <- df_rcll_l_sub %>% 
  rstatix::t_test(s_RT ~ suppression, paired = F, var.equal = F)

l_sub_rcll_l <- df_rcll_l_sub %>% 
  tidyr::drop_na() %>% 
  split(.$suppression)
g_rcll_l_sub <- MBESS::smd(l_sub_rcll_l$`0`$s_RT, l_sub_rcll_l$`1`$s_RT, Unbiased = T)

# reporting
format_aov(aov_rcll_d_e2)
format_aov(sme_rcll_d_e2, grouped = T)
dplyr::left_join(mc_rcll_d_e2, df_es_rcll_d_e2, by = c("level1", "level2")) %>% 
  format_mc(es = T)
ttest_sub_rcll_l_e2["g"] <- g_rcll_l_sub #, df_es_nt_ques_e1, by = c("suppression", "group1" = "level1", "group2" = "level2")) %>% 
format_t(ttest_sub_rcll_l_e2, es = T)

# dot probe ---------------------------------------------------------------

df_dp_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, matches("(cong|incong)$"), new) %>% 
  tidyr::pivot_longer(matches("(cong|incong)$"), 
                      names_to = c("condition", "type"),
                      names_sep = "_",
                      values_to = "latency") %>% 
  dplyr::mutate(diff = (latency - new) / new)

# anova
model_dp_e2 <- diff ~ suppression*condition*type + Error(participant/(condition*type))
aov_dp_e2 <- df_dp_e2 %>% 
  dplyr::select(-new, -latency) %>% 
  rstatix::anova_test(model_dp_e2) %>% 
  rstatix::get_anova_table(correction = "auto")

# simple main effect
sme_dp_e2 <- df_dp_e2 %>% 
  dplyr::select(-new, -latency) %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(diff ~ condition*type + Error(participant/(condition*type))) %>% 
  rstatix::get_anova_table(correction = "auto")

# reporting
format_aov(aov_dp_e2)
format_aov(sme_dp_e2, grouped = T)

# correlation -------------------------------------------------------------

df_corr <- df_e2 %>% 
  dplyr::select(participant, suppression, ends_with("cong"), ends_with("recall"), ends_with("delay")) %>% 
  tidyr::pivot_longer(cols = c(-participant, -suppression),
                      names_to = c("TNTcond", "variable"),
                      names_sep = "_",
                      values_to = "value") %>% 
  dplyr::filter(TNTcond != "s") %>% 
  tidyr::pivot_wider(names_from = variable, values_from = value)

calc_scor <- function(df, group){
  df %>% 
    dplyr::filter(suppression == group) %>% 
    tidyr::drop_na() %>% 
    split(.$TNTcond) %>% 
    purrr::map(.f = dplyr::select, delay, cong, incong) %>% 
    purrr::map(mscorci)
}

res_corr_ds <- calc_scor(df_corr, group = 0)
res_corr_ts <- calc_scor(df_corr, group = 1)



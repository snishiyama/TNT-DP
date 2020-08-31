
# setup -------------------------------------------------------------------

if (!exists("df_e1")) {
  source(here::here("analyses/setup.R"))
}


# questionnaire -----------------------------------------------------------

df_ques_e1 <- df_e1 %>% 
  dplyr::select(participant, suppression, starts_with("think"), starts_with("nothink")) %>% 
  tidyr::pivot_longer(cols = -c(suppression, participant), 
                      names_to = c("condition", "type"), 
                      names_sep = "_",
                      values_to = "rating")
# ANOVA
model_ques_e1 <- rating ~ suppression*type+Error(participant/type)
l_aov_ques_e1 <- df_ques_e1 %>% 
  split(.$condition) %>% 
  purrr::map(.f = rstatix::anova_test, model_ques_e1) %>% 
  purrr::map(.f = rstatix::get_anova_table, correction = "auto")

# multiple comparisons for cue types of think items
mc_ques_think_e1 <- multi_aov(l_aov_ques_e1$think, "type")

# group-by simple main effect analyses for no-think items
sme_ques_nt_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(rating ~ type + Error(participant/type)) %>% 
  rstatix::get_anova_table(correction = "auto")

# pairwise t test 
pwt_ques_nt_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)
df_es_nt_ques_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "nothink") %>% 
  tidyr::pivot_wider(names_from = type, values_from = rating) %>% 
  group_by(suppression) %>% 
  dplyr::summarise(cue_sub = MBESS::smd(cue, sub, Unbiased = T),
                   cue_tar = MBESS::smd(cue, tar, Unbiased = T),
                   sub_tar = MBESS::smd(tar, sub, Unbiased = T)) %>% 
  tidyr::pivot_longer(cols = cue_sub:sub_tar, names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# reporting
purrr::map(l_aov_ques_e1, format_aov)
format_mc(mc_ques_think_e1, es = T)
dplyr::left_join(pwt_ques_nt_e1, df_es_nt_ques_e1, by = c("suppression", "group1" = "level1", "group2" = "level2")) %>% 
  format_t(grouped = T, es = T)

# recall rate -------------------------------------------------------------

df_rcll_e1 <- df_e1 %>% 
  dplyr::select(participant, suppression, 
                think = recall_think, 
                nothink = recall_nothink, 
                baseline = recall_base, 
                substitute = recall_sub) %>% 
  tidyr::pivot_longer(think:substitute, names_to = "status", values_to = "rate")

# anova t vs nt vs b
model_rcll_e1 <- rate ~ suppression*status + Error(participant/status)
aov_rcll_e1 <- df_rcll_e1 %>% 
  dplyr::filter(status != "substitute") %>% 
  rstatix::anova_test(model_rcll_e1) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_e1 <- multi_aov(aov_rcll_e1, "status")

# source(here::here("analyses/anovakun_485.txt"))
# aovkun <- df_rcll_e1 %>%
#   dplyr::filter(status != "substitute") %>%
#   anovakun("AsB", long = T, geta = T)

# multiple comparisons with emmeans
# df_rcll_e1 %>% 
#   dplyr::filter(status != "substitute") %>% 
#   dplyr::mutate(participant = as_factor(participant),
#                 suppression = as_factor(suppression),
#                 status = factor(status, levels = c("think", "nothink", "baseline"))) %>% 
#   aov(rate ~ suppression*status+Error(participant/status), 
#       contrasts = list(suppression = "contr.sum", status = "contr.sum"), # type 3 SS
#       data = .) %>% 
#   emmeans::emmeans(~ status) %>% 
#   emmeans::contrast(method = "pairwise", adjust = "holm")

# t.test on substitute
ttest_sub_rcll_e1 <- df_rcll_e1 %>% 
  dplyr::filter(status == "substitute") %>% 
  rstatix::t_test(rate ~ suppression)

l_sub_rcll_e1 <- df_rcll_e1 %>% 
    dplyr::filter(status == "substitute") %>% 
    split(.$suppression)

ttest_sub_rcll_e1["g"] <- MBESS::smd(l_sub_rcll_e1$`0`$rate, l_sub_rcll_e1$`1`$rate)

# reporting
format_aov(aov_rcll_e1)
format_mc(mc_rcll_e1, es = T)
format_t(ttest_sub_rcll_e1, es = T, p_adj = F)

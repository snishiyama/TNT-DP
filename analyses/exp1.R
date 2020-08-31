
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
mc_ques_think_e1 <- df_ques_e1 %>% 
  dplyr::filter(condition == "think") %>% 
  multi_aov("type", model_ques_e1, c("cue", "tar", "sub"))

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
mc_rcll_e1 <- multi_aov(df_rcll_e1, "status", model_rcll_e1, c("think", "nothink", "baseline"))

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

# Effect size
df_es_rcll_e1 <- df_rcll_e1 %>% 
  tidyr::pivot_wider(names_from = status, values_from = rate) %>% 
  dplyr::summarise(g_tnt = MBESS::smd(think, nothink, Unbiased = T),
                   g_tb = MBESS::smd(think, baseline, Unbiased = T),
                   g_ntb = MBESS::smd(nothink, baseline, Unbiased = T))

# t.test on substitute
ttest_sub_rcll_e1 <- df_rcll_e1 %>% 
  dplyr::filter(status == "substitute") %>% 
  rstatix::t_test(rate ~ suppression)

l_sub_rcll_e1 <- df_rcll_e1 %>% 
    dplyr::filter(status == "substitute") %>% 
    split(.$suppression)

g_rcll_sub_e1 <- MBESS::smd(l_sub_rcll_e1$`0`$rate, l_sub_rcll_e1$`1`$rate)



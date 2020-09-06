
# setup -------------------------------------------------------------------

if (!exists("df_e1")) {
  source(here::here("analyses/setup.R"))
}


# questionnaire -----------------------------------------------------------

df_ques_mean_e1 <- df_ques_e1 %>% 
  dplyr::group_by(participant, suppression, status, type) %>% 
  dplyr::summarise(rating = mean(rate), .groups = "drop") %>% 
  dplyr::mutate_if(is.character, as_factor)

# Table
# df_ques_mean_e1 %>% 
#   dplyr::mutate(status = fct_relevel(status, "think"),
#                 type = fct_relevel(type, "cue", "tar")) %>% 
#   dplyr::group_by(suppression, status, type) %>% 
#   dplyr::summarise(value = mean_cl(rating), name = c("Mean", "CI"), .groups = "drop") %>% 
#   tidyr::unite(col = "cond", suppression, type, sep = "_") %>% 
#   tidyr::pivot_wider(names_from = cond, values_from = value)

# ANOVA
l_aov_ques_e1 <- df_ques_mean_e1 %>% 
  split(.$status) %>% 
  purrr::map(.f = rstatix::anova_test, rating ~ suppression*type+Error(participant/type)) %>% 
  purrr::map(.f = rstatix::get_anova_table, correction = "auto")

# multiple comparisons for cue types of think items
mc_ques_think_e1 <- multi_aov(l_aov_ques_e1$think, "type")

# group-by simple main effect analyses for no-think items
sme_ques_nt_e1 <- df_ques_mean_e1 %>% 
  dplyr::filter(status == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(rating ~ type + Error(participant/type)) %>% 
  rstatix::get_anova_table(correction = "auto")

# pairwise t test 
pwt_ques_nt_e1 <- df_ques_mean_e1 %>% 
  dplyr::filter(status == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)
df_es_nt_ques_e1 <- df_ques_mean_e1 %>% 
  dplyr::filter(status == "nothink") %>% 
  tidyr::pivot_wider(names_from = type, values_from = rating) %>% 
  dplyr::group_by(suppression) %>% 
  dplyr::summarise(cue_sub = MBESS::smd(cue, sub, Unbiased = T),
                   cue_tar = MBESS::smd(cue, tar, Unbiased = T),
                   sub_tar = MBESS::smd(tar, sub, Unbiased = T),
                   .groups = "drop") %>% 
  tidyr::pivot_longer(cols = -suppression, names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# reporting
purrr::map(l_aov_ques_e1, format_aov)
format_mc(mc_ques_think_e1, es = T)
format_aov(sme_ques_nt_e1, grouped = T)
dplyr::left_join(pwt_ques_nt_e1, df_es_nt_ques_e1, by = c("suppression", "group1" = "level1", "group2" = "level2")) %>% 
  format_t(grouped = T, es = T)


# recall rate -------------------------------------------------------------

# df_rcll_e1 <- df_e1 %>% 
#   dplyr::select(participant, suppression, 
#                 think = recall_think, 
#                 nothink = recall_nothink, 
#                 baseline = recall_base, 
#                 substitute = recall_sub) %>% 
#   tidyr::pivot_longer(think:substitute, names_to = "status", values_to = "rate")

df_rcll_e1 <- df_e1 %>% 
  dplyr::select(participant, suppression, item_id = ID, status = condition,
                pre_recall = learning, post_recall = target_recall) %>% 
  dplyr::filter(pre_recall == 1) %>% 
  dplyr::group_by(participant, suppression, status) %>% 
  dplyr::summarise(rate = mean(post_recall), .groups = "drop")

# anova t vs nt vs b
aov_rcll_e1 <- df_rcll_e1 %>% 
  rstatix::anova_test(rate ~ suppression*status + Error(participant/status)) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_e1 <- multi_aov(aov_rcll_e1, "status")

# t.test on substitute
df_rcll_sub_e1 <- df_e1 %>% 
  dplyr::select(participant, ID, suppression, sub_recall) %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(participant, suppression) %>% 
  dplyr::summarise(rate = mean(sub_recall), .groups = "drop")

ttest_sub_rcll_e1 <- df_rcll_sub_e1 %>% 
  rstatix::t_test(rate ~ suppression)

l_sub_rcll_e1 <- df_rcll_sub_e1 %>%
    split(.$suppression)

ttest_sub_rcll_e1["g"] <- MBESS::smd(l_sub_rcll_e1$DS$rate, l_sub_rcll_e1$TS$rate, Unbiased = T)

# reporting
format_aov(aov_rcll_e1)
format_mc(mc_rcll_e1, es = T)
format_t(ttest_sub_rcll_e1, es = T, p_adj = F)

# Table
# df_rcll_e1 %>%
#   dplyr::bind_rows(df_rcll_sub_e1 %>% dplyr::mutate(status = "substitute")) %>% 
#   dplyr::mutate(status = factor(status, levels = c("think", "nothink", "baseline", "substitute")),
#                 suppression = factor(suppression, levels = c("DS", "TS"))) %>% 
#   dplyr::group_by(suppression, status) %>%
#   dplyr::summarise(value = mean_cl(rate), name = c("Mean", "CI"), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = status, values_from = value)

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
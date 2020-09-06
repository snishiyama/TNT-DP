
# setup -------------------------------------------------------------------

if (!exists("df_e2")) {
  source(here::here("analyses/setup.R"))
}


# questionnaire -----------------------------------------------------------

df_ques_mean_e2 <- df_ques_e2 %>% 
  dplyr::group_by(participant, suppression, status, type) %>% 
  dplyr::summarise(rating = mean(rate), .groups = "drop") %>% 
  dplyr::mutate_if(is.character, as_factor)

# ANOVA
l_aov_ques_e2 <- df_ques_mean_e2 %>% 
  split(.$status) %>% 
  purrr::map(.f = rstatix::anova_test, rating ~ suppression*type+Error(participant/type)) %>% 
  purrr::map(.f = rstatix::get_anova_table, correction = "auto")

# multiple comparisons for cue types of think items
mc_ques_think_e2 <- multi_aov(l_aov_ques_e2$think, "type")

# group-by simple main effect analyses for no-think items
sme_ques_nt_e2 <- df_ques_mean_e2 %>% 
  dplyr::filter(status == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(rating ~ type + Error(participant/type)) %>% 
  rstatix::get_anova_table(correction = "auto")

# pairwise t test 
pwt_ques_nt_e2 <- df_ques_mean_e2 %>% 
  dplyr::filter(status == "nothink") %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::pairwise_t_test(rating ~ type, paired = T)

# effect sizes
df_es_nt_ques_e2 <- df_ques_mean_e2 %>% 
  dplyr::filter(status == "nothink") %>% 
  tidyr::pivot_wider(names_from = type, values_from = rating) %>% 
  dplyr::group_by(suppression) %>% 
  dplyr::summarise(cue_sub = MBESS::smd(cue, sub, Unbiased = T),
                   cue_tar = MBESS::smd(cue, tar, Unbiased = T),
                   sub_tar = MBESS::smd(tar, sub, Unbiased = T), .groups = "drop") %>% 
  tidyr::pivot_longer(cols = -suppression, names_to = c("level1", "level2"), names_sep = "_", values_to = "g")

# reporting
purrr::map(l_aov_ques_e2, format_aov)
format_mc(mc_ques_think_e2, es = T)
format_aov(sme_ques_nt_e2, group = T)
dplyr::left_join(pwt_ques_nt_e2, df_es_nt_ques_e2, by = c("suppression", "group1" = "level1", "group2" = "level2")) %>% 
  format_t(grouped = T, es = T)


# recall rate -------------------------------------------------------------

df_rcll_r_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, item_id, status, contains("recall")) %>% #post_recall_corr) %>%
  dplyr::filter(pre_recall_corr == 1) %>% 
  dplyr::group_by(participant, suppression, status) %>% 
  dplyr::summarise(rate = mean(post_recall_corr), .groups = "drop")

# anova t vs nt vs b
aov_rcll_r_e2 <- df_rcll_r_e2 %>% 
  rstatix::anova_test(rate ~ suppression*status + Error(participant/status)) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_r_e2 <- multi_aov(aov_rcll_r_e2, "status")

# t.test on substitute
df_sub_rcll_r_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, item_id, sub_recall) %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(participant, suppression) %>% 
  dplyr::summarise(rate = mean(sub_recall), .groups = "drop")

ttest_sub_rcll_r_e2 <- df_sub_rcll_r_e2 %>% 
  rstatix::t_test(rate ~ suppression, paired = F, var.equal = F)

l_sub_rcll_r_e2 <- df_sub_rcll_r_e2 %>% 
  split(.$suppression)

ttest_sub_rcll_r_e2["g"] <- MBESS::smd(l_sub_rcll_r_e2$DS$rate, l_sub_rcll_r_e2$TS$rate, Unbiased = T)

# reporting
format_aov(aov_rcll_r_e2)
format_mc(mc_rcll_r_e2, es = T)
format_t(ttest_sub_rcll_r_e2, p_adj = F, es = T)

# Table
df_rcll_r_e2 %>%
  dplyr::bind_rows(df_sub_rcll_r_e2 %>% dplyr::mutate(status = "substitute")) %>%
  dplyr::mutate(status = factor(status, levels = c("think", "nothink", "baseline", "substitute")),
                suppression = factor(suppression, levels = c("DS", "TS"))) %>%
  dplyr::group_by(suppression, status) %>%
  dplyr::summarise(value = mean_cl(rate), name = c("Mean", "CI"), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = status, values_from = value)

# recall latency ----------------------------------------------------------

# df_rcll_d_e2 <- df_e2 %>% 
#   dplyr::select(participant, suppression, ends_with("delay")) %>% 
#   dplyr::rename(think = t_delay, nothink = nt_delay, baseline = b_delay) %>% 
#   tidyr::pivot_longer(think:baseline, names_to = "status", values_to = "delay") %>% 
#   tidyr::drop_na() # remove data of id 27

df_rcll_d_e2 <- df_e2 %>% 
  dplyr::select(participant, suppression, item_id, status, contains("recall")) %>% 
  dplyr::filter(pre_recall_corr == 1, post_recall_corr == 1) %>% 
  dplyr::mutate(diff = post_recall_rt - pre_recall_rt) %>% 
  dplyr::group_by(participant, suppression, status) %>% 
  dplyr::summarise(delay = mean(diff, na.rm = T), .groups = "drop")

# anova t vs nt vs b
aov_rcll_d_e2 <- df_rcll_d_e2 %>% 
  rstatix::anova_test(delay ~ suppression*status + Error(participant/status)) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_d_e2 <- multi_aov(aov_rcll_d_e2, "status")

# group-by simple main effect analyses for no-think items
sme_rcll_d_e2 <- df_rcll_d_e2 %>% 
  dplyr::group_by(status) %>% 
  rstatix::anova_test(delay ~ suppression) %>% 
  rstatix::get_anova_table(correction = "auto")

# substitute
df_rcll_l_sub_mean <- df_rcll_l_sub %>% 
  dplyr::filter(sub_corr != 0) %>% 
  dplyr::group_by(participant, suppression) %>% 
  dplyr::summarise(mean = mean(sub_rt, na.rm = T), .groups = "drop")

ttest_sub_rcll_l_e2 <- df_rcll_l_sub_mean %>% 
  rstatix::t_test(mean ~ suppression, paired = F, var.equal = F)

l_sub_rcll_l <- df_rcll_l_sub_mean %>% 
  split(.$suppression)

ttest_sub_rcll_l_e2["g"] <- MBESS::smd(l_sub_rcll_l$DS$mean, l_sub_rcll_l$TS$mean, Unbiased = T)

# reporting
format_aov(aov_rcll_d_e2)
format_mc(mc_rcll_d_e2, es = T)
format_aov(sme_rcll_d_e2, grouped = T)
format_t(ttest_sub_rcll_l_e2, p_adj = F, es = T)


# remove id 15, 27 ---

df_rcll_d_e2_rm <- df_rcll_d_e2 %>% dplyr::filter(!participant %in% c(15, 27))

# anova t vs nt vs b
aov_rcll_d_e2_rm <- df_rcll_d_e2_rm %>% 
  rstatix::anova_test(delay ~ suppression*status + Error(participant/status)) %>% 
  rstatix::get_anova_table(correction = "auto")

# multiple comparisons
mc_rcll_d_e2_rm <- multi_aov(aov_rcll_d_e2_rm, "status")

# group-by simple main effect analyses for no-think items
sme_rcll_d_e2_rm <- df_rcll_d_e2_rm %>% 
  dplyr::group_by(status) %>% 
  rstatix::anova_test(delay ~ suppression) %>% 
  rstatix::get_anova_table(correction = "auto")

format_aov(aov_rcll_d_e2_rm)
format_mc(mc_rcll_d_e2_rm, es = T)
format_aov(sme_rcll_d_e2_rm, grouped = T)


# dot probe ---------------------------------------------------------------

df_dp_e2 <- df_dp_raw_e2 %>% 
  filter(
    ((participant < 25 & corr == 1) | (participant >= 25 & corr == 0))#,
    # RT < 1
  ) %>% 
  dplyr::group_by(participant) %>% 
  dplyr::mutate(mad = mad(RT),
                median = median(RT),
                suppression = if_else(participant %% 2 == 0, "TS", "DS")) %>% 
  dplyr::filter(between(RT, median - 2.5 * mad, median + 2.5 * mad))

df_dp_rm <- dplyr::anti_join(df_dp_raw_e2, 
                 df_dp_e2 %>% dplyr::select(participant, ID), 
                 by = c("participant", "ID")) %>% 
  dplyr::mutate(filler = if_else(status == "filler", 1, 0)) %>% 
  dplyr::group_by(participant, filler) %>% 
  dplyr::summarise(freq = n(), .groups = "drop") %>% 
  dplyr::group_by(filler) %>%
  dplyr::summarise(mean = mean(freq), sd = sd(freq), .groups = "drop")

df_dp_rm$mean / c(36, 124)

df_dp_mean_e2 <- df_dp_e2 %>% 
  dplyr::left_join(df_e2 %>% dplyr::select(participant, item_id, status, pre_recall_corr),
                   by = c("participant", "ID" = "item_id", "status")) %>% 
  dplyr::filter(pre_recall_corr == 1) %>% 
  dplyr::group_by(participant, suppression, status, congruency) %>% 
  dplyr::summarise(mean = mean(RT), .groups = "drop") %>% 
  dplyr::filter(status != "filler")

# anova
aov_dp_e2 <- df_dp_mean_e2 %>% 
  rstatix::anova_test(mean ~ suppression*status*congruency + Error(participant/(status*congruency))) %>% 
  rstatix::get_anova_table(correction = "auto")

# simple main effect
sme_dp_e2 <- df_dp_mean_e2 %>% 
  dplyr::group_by(suppression) %>% 
  rstatix::anova_test(mean ~ status*congruency + Error(participant/(status*congruency))) %>% 
  rstatix::get_anova_table(correction = "auto")

# reporting
format_aov(aov_dp_e2)
format_aov(sme_dp_e2, grouped = T)

# trial-by-trial variances
df_dp_e2 %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(participant, sd) %>% 
  dplyr::summarise(mean = mean(sd))

# correlation -------------------------------------------------------------

df_corr <- dplyr::left_join(
    df_rcll_d_e2,
    tidyr::pivot_wider(df_dp_mean_e2, names_from = "congruency", values_from = "mean"),
    by = c("participant", "suppression", "status")
  ) %>% 
  dplyr::mutate(dp = incong - cong)

calc_scor <- function(df, group) {
  df %>% 
    dplyr::filter(suppression == group) %>% 
    tidyr::drop_na() %>% 
    split(.$status) %>% 
    purrr::map(.f = dplyr::select, delay, cong, incong, dp) %>% 
    purrr::map(mscorci, corfun = pcor, nboot = 500)
}

# res_corr_ds <- calc_scor(df_corr, group = "DS")
res_corr_ts <- calc_scor(df_corr, group = "TS")

df_corr_rm <- df_corr %>% dplyr::filter(!participant %in% c(15, 27))

res_corr_ds_rm <- calc_scor(df_corr_rm, group = "DS")

df_corr_rm %>%
  dplyr::filter(suppression == "DS") %>%
  split(.$status) %>% 
  purrr::map(dplyr::select, delay, cong, incong, dp) %>%
  purrr::map(psych::corr.test, method = "pearson")

df_corr_rm %>%
  dplyr::filter(suppression == "TS") %>%
  split(.$status) %>% 
  purrr::map(dplyr::select, delay, cong, incong, dp) %>%
  purrr::map(psych::corr.test, method = "pearson")

# 
df_corr_rm %>%
  dplyr::filter(suppression == "DS", status == "nothink") %>%
  ggplot() +
  aes(x = delay, y = dp*1000) +
  geom_point()

# by-item correlation --
df_dp_delay <- dplyr::left_join(
  x = df_e2 %>% 
    dplyr::filter(pre_recall_corr == 1, post_recall_corr == 1) %>% 
    dplyr::mutate(diff = post_recall_rt - pre_recall_rt) %>% 
    dplyr::select(participant, suppression, item_id, status, diff),
  y = df_dp_e2 %>% dplyr::select(participant, suppression, status, congruency, "item_id" = ID, "dp_rt" = RT),
  by = c("participant", "suppression", "status", "item_id"))

df_dp_delay %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(participant, suppression, status, congruency) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(corr = purrr::map(.x = data, .f = ~cor(.x$diff, .x$dp_rt, method = "spearman"))) %>% 
  tidyr::unnest(corr) %>%
  dplyr::group_by(suppression, status, congruency) %>% 
  dplyr::summarise(mean = mean(corr, na.rm = T), sd = sd(corr, na.rm = T))
  

# plot --------------------------------------------------------------------

gg_DP <- df_dp_mean_e2 %>% 
  mutate(suppression = if_else(suppression == "DS", "Direct suppression", "Thought substitution"),
         type = if_else(congruency == "cong", "Congruent", "Incongruent"),
         condition = case_when(status == "think" ~ "Think",
                               status == "nothink" ~ "No-Think",
                               status == "baseline" ~ "Baseline"),
         condition = forcats::fct_relevel(condition, "Think", "No-Think")) %>% 
  ggplot(aes(x = condition, y = mean, fill = type)) +
  stat_summary(geom = "bar", fun = mean, position = position_dodge2(), color = "black", size = 0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal, 
               position = position_dodge(width = 0.9), size = 0.3, width = 0.3) +
  scale_fill_grey(start = 1, end = 0.6) +
  labs(y = "RT (s)") +
  # geom_hline(yintercept = 0, size = 0.3) +
  coord_cartesian(ylim = c(0.45, 0.55)) +
  facet_wrap(~ suppression, ncol = 2, strip.position = "bottom") +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(margin = margin(t = 0, r = 8, b = 0, l = 3, unit = "pt")),
        legend.position = "top", 
        strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text = element_text(size=12)
  )

gg_delay <- df_rcll_d_e2 %>% 
  mutate(suppression = if_else(suppression == "DS", "Direct suppression", "Thought substitution"),
         condition = case_when(status == "think" ~ "Think",
                               status == "nothink" ~ "No-Think",
                               status == "baseline" ~ "Baseline"),
         condition = forcats::fct_relevel(condition, "Think", "No-Think")) %>% 
  ggplot(aes(x = suppression, y = delay, fill = condition)) +
  stat_summary(geom = "bar", fun = mean, position = position_dodge2(), color = "black", size = 0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal, 
               position = position_dodge(width = 0.9), size = 0.3, width = 0.3) +
  scale_fill_grey(start = 1, end = 0.6) +
  geom_hline(yintercept = 0, size = 0.3) +
  labs(y = "Recall delay (s)") +
  # facet_grid(cols = vars(suppression)) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(t = 0, r = 8, b = 0, l = 2, unit = "pt")),
        legend.position = "top", 
        legend.background = element_blank())

# ggsave(filename = here::here("analyses/figure/gg_DP_exp2.pdf"), plot = gg_DP, height = 90, width = 160, units = "mm", dpi = 300)
# ggsave(filename = here::here("analyses/figure/gg_delay_exp2.pdf"), plot = gg_delay, height = 90, width = 160, units = "mm", dpi = 300)

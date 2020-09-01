
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
  dplyr::group_by(participant, suppression) %>% 
  dplyr::summarise(mean = mean(rt, na.rm = T), .groups = "drop")

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


# dot probe ---------------------------------------------------------------

readr::read_csv(here::here("exp/data/TNT_exp2_sotsuron.csv")) %>%
  dplyr::mutate(suppression = suppression - 1L) %>%
  dplyr::select(participant, suppression, matches("(cong|incong)$"), new) %>% 
  tidyr::pivot_longer(matches("(cong|incong)$"), 
                      names_to = c("condition", "type"),
                      names_sep = "_",
                      values_to = "latency") %>% 
  ggplot() +
  aes(x = condition, y = latency, color = type) +
  stat_summary(geom = "pointrange", fun.data = "mean_se", position = position_dodge(width = 0.5)) +
  facet_wrap(~suppression)
# %>% 
#   dplyr::mutate(diff = (latency - new) / new)
readr::read_csv(here::here("exp/data/TNT_exp2_sotsuron.csv")) %>%
  dplyr::mutate(suppression = suppression - 1L) %>%
  dplyr::select(participant, suppression, matches("(cong|incong)$"), new) %>% 
  tidyr::pivot_longer(matches("(cong|incong)$"), 
                      names_to = c("condition", "type"),
                      names_sep = "_",
                      values_to = "latency") %>% 
  dplyr::select(-new) %>%
  rstatix::anova_test(latency ~ suppression*condition*type + Error(participant/(condition*type))) %>% 
  rstatix::get_anova_table(correction = "auto")

df_dp_e2 <- df_e2 %>% 
  dplyr::filter(pre_recall_corr == 1, 
                ((participant < 25 & dp_corr == 1) | (participant >= 25 & dp_corr == 0)),
                dp_rt < 1
                ) %>% 
  dplyr::group_by(participant) %>% 
  dplyr::mutate(mean = mean(dp_rt), sd = sd(dp_rt)) %>% 
  dplyr::filter(between(dp_rt, mean - 2.5 * sd, mean + 2.5 * sd))

df_dp_mean_e2 <- df_dp_e2 %>% 
  dplyr::group_by(participant, suppression, status, congruency) %>% 
  dplyr::summarise(mean = mean(dp_rt), .groups = "drop")

df_dp_mean_e2 %>% 
  ggplot() +
  aes(x = status, y = mean, color = congruency) +
  stat_summary(geom = "pointrange", fun.data = "mean_se", position = position_dodge(width = 0.5)) +
  facet_wrap(~suppression)

# anova
# model_dp_e2 <- diff ~ suppression*condition*type + Error(participant/(condition*type))
model_dp_e2 <- mean ~ suppression*status*congruency + Error(participant/(status*congruency))
aov_dp_e2 <- df_dp_mean_e2 %>% 
  # dplyr::select(-new, -latency) %>% 
  rstatix::anova_test(mean ~ suppression*status*congruency + Error(participant/(status*congruency))) %>% 
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

calc_scor <- function(df, group) {
  df %>% 
    dplyr::filter(suppression == group) %>% 
    tidyr::drop_na() %>% 
    split(.$TNTcond) %>% 
    purrr::map(.f = dplyr::select, delay, cong, incong) %>% 
    purrr::map(mscorci)
}

res_corr_ds <- calc_scor(df_corr, group = 0)
res_corr_ts <- calc_scor(df_corr, group = 1)


# plot --------------------------------------------------------------------

gg_DP <- df_dp_e2 %>% 
  mutate(suppression = if_else(suppression == 0, "Direct suppression", "Thought substitution"),
         type = if_else(type == "cong", "Congruent", "Incongruent"),
         condition = case_when(condition == "t" ~ "Think",
                               condition == "nt" ~ "No-Think",
                               condition == "b" ~ "Baseline"),
         condition = forcats::fct_relevel(condition, "Think", "No-Think")) %>% 
  ggplot(aes(x = condition, y = diff, fill = type)) +
  stat_summary(geom = "bar", fun = mean, position = position_dodge2(), color = "black", size = 0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal, 
               position = position_dodge(width = 0.9), size = 0.3, width = 0.3) +
  scale_fill_grey(start = 1, end = 0.6) +
  labs(y = "Adjusted RT scores") +
  geom_hline(yintercept = 0, size = 0.3) +
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
  labs(y = "Recall delay (ms)") +
  # facet_grid(cols = vars(suppression)) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(t = 0, r = 8, b = 0, l = 2, unit = "pt")),
        legend.position = "top", 
        legend.background = element_blank())


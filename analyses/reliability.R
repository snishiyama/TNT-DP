
# setup -------------------------------------------------------------------

if (!exists("df_e1")) {
  source(here::here("analyses/setup.R"))
}

# check onset reliability -------------------------------------------------

get_sheet_regexp <- function(filepath, pattern, cols = NULL) {
  par_id <- stringr::str_extract(filepath, "\\d+_TNT") %>% readr::parse_number()
  # print(par_id)
  snm <- readxl::excel_sheets(filepath) %>% str_subset(pattern = pattern)
  df <- NULL
  tryCatch({
    df <- readxl::read_xlsx(path = filepath, sheet = snm, col_names = TRUE) %>%
      dplyr::mutate(participant = par_id)
    if (!is.null(cols)) {
      df <- df %>% dplyr::select(participant, !!cols) %>% tidyr::drop_na()
    }
    return(df)
  },
  error = function(e){
    message(e)
  })
}

df_pre_e2 <- purrr::map_dfr(
  .x = fs::dir_ls(here::here("exp/data/exp2"), glob = "*.xlsx"),
  .f = get_sheet_regexp, pattern = "learn.*(test|Test)", cols = c("ID", "order")
) %>% dplyr::mutate(phase = "pre", order = order + 1)

df_post_e2 <- purrr::map_dfr(
  .x = fs::dir_ls(here::here("exp/data/exp2"), glob = "*.xlsx"),
  .f = get_sheet_regexp, pattern = "tar.*recall", cols = c("ID", "order")
) %>% dplyr::mutate(phase = "post", order = order + 1)

df_sub_e2 <- purrr::map_dfr(
  .x = fs::dir_ls(here::here("exp/data/exp2"), glob = "*.xlsx"),
  .f = get_sheet_regexp, pattern = "sub_(r|R)ecall", cols = c("ID", "order")
) %>% dplyr::mutate(phase = "sub", order = order + 1)

df_chronset <- fs::dir_ls(here::here("exp/data/chronset")) %>% 
  purrr::map_dfr(function(x){
    par_id <- stringr::str_extract(x, "\\d+_TNT") %>% readr::parse_number()
    df <- readr::read_tsv(x, col_names = c("filename", "chronset")) %>%
      dplyr::mutate(participant = par_id)
    return(df)
  }) %>% 
  dplyr::mutate(phase = case_when(str_detect(filename, "learn") ~ "pre",
                                  str_detect(filename, "sub") ~ "sub",
                                  TRUE ~ "post"),
                order = readr::parse_number(filename)) %>% 
  dplyr::left_join(bind_rows(df_pre_e2, df_post_e2, df_sub_e2), 
                   by = c("participant", "phase", "order"))


df_manual <- df_e2 %>% 
  dplyr::select(participant, "ID" = item_id, pre_recall_rt, post_recall_rt) %>% 
  tidyr::pivot_longer(cols = contains("recall_rt"), names_to = "phase", values_to = "manual") %>% 
  dplyr::mutate(phase = str_extract(phase, "(pre|post)")) %>% 
  dplyr::bind_rows(df_rcll_l_sub %>% 
                     dplyr::select(participant, ID, "manual" = sub_rt) %>%
                     dplyr::mutate(phase = "sub")) %>% 
  dplyr::mutate(manual = manual * 1000)

df_chr_man <- dplyr::left_join(df_chronset, df_manual, by = c("participant", "phase", "ID")) %>% 
  tidyr::drop_na()

df_record_w_filler <- readxl::read_xlsx(here::here("exp/data/record_w_filler.xlsx")) %>% 
  dplyr::select(participant, phase, order, ID, filler)

df_chr_man_new <- left_join(df_chr_man, df_record_w_filler, by = c("participant", "phase", "order", "ID")) %>%
  tidyr::replace_na(list(filler = 0)) %>% 
  dplyr::filter(filler != 1) %>% 
  dplyr::filter(chronset > 500)

res_lm <- lm(chronset ~ manual, data = df_chr_man_new)
R2 <- (var(res_lm$fitted.values) * res_lm$df.residual) / (var(res_lm$model$chronset) * res_lm$df.residual)
gg_chr_man <- df_chr_man_new %>% 
  ggplot2::ggplot() +
  aes(x = manual, y = chronset) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  geom_point(shape = 1, size = 1.5, alpha = 0.5) +
  stat_smooth(formula = y ~ x, method = "lm", color = "black", fullrange = T, se = F) +
  coord_cartesian(xlim = c(0, 3200), ylim = c(0, 3200), expand = F) +
  xlim(0,3200) +
  annotate("text", 
           label = sprintf("paste(y == %.2f + %.2f * x, ','~~italic(R) ^ 2==%.2f)", res_lm$coefficients["(Intercept)"], res_lm$coefficients["manual"], R2), 
           x = 1200, y = 2800, parse = T) +
  theme_bw(base_family = "Helvetica") +
  labs(x = "Manual (ms)",y = "Chronset (ms)") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))

# cor(df_chr_man_new$manual, df_chr_man_new$chronset)
# ggsave(here::here("analyses/figure/gg_chr_man.pdf"), gg_chr_man, height = 90, width = 100, units = "mm", dpi = 300)
# df_chr_man %>%
#   dplyr::filter(manual - 300 > chronset, chronset > 500) %>%
#   dplyr::arrange(participant) %>%
#   write_csv(here::here("exp/data/df_man-over-chr.csv"))
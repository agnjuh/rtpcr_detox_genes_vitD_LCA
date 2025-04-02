library(magrittr)
library(readxl)
library(dplyr)
library(rlang)
library(stringr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)

CONTROLS <- c('ET', 'MET')
MET_CONTROL <- c('RIF')
MET_CONTROL %<>% c('MET')
DATASETS <- c('1', '2')
SINGLE_TBR_PLOT <- c('CYP2B6')
TS_RANGES <- c(h12 = 'A1:H64', h24 = 'K1:R65', h48 = 'T1:AA65')
XLABELS_1_2 <- c('ET', '1,25D', '3KLCA', 'MET', 'RIF')
XLABELS_3 <- c('ET', '1,25D', 'EB1089', '3KLCA', 'MET', 'RIF')

rtpcr_main <- function() {
  
  s <-
    DATASETS %>%
    purrr::map(
      rtpcr_workflow
    ) %T>%
    {rlang::exec(rtpcr_summary, !!!.)}
  
  t <-
    read_time_series('3') %T>%
    rtpcr_export() %T>%
    rtpcr_export_nice()
  
  rtpcr_times_plot(t, fold_induction, ylabel = 'Fold induction')
  rtpcr_times_plot(t, targ_by_ref_d, ylabel = 'Target to reference ratio')
  
  rtpcr_summary_t(t)
  
  invisible(list(single = s, times = t))
  
}


rtpcr_workflow <- function(dataset) {
  
  d <- read_and_preprocess(dataset)
  
  rtpcr_barplot(d, targ_by_ref_d, ylabel = 'Target to reference ratio')
  rtpcr_barplot(d, fold_induction, ylabel = 'Fold induction')
  
  d %>%
    dplyr::pull(target) %>%
    intersect(SINGLE_TBR_PLOT) %>%
    purrr::walk(
      function(protein) {
        d %>%
          filter(target == protein) %>%
          `attr<-`('idx', protein) %>%
          rtpcr_barplot(targ_by_ref_d, ylabel = 'Target to reference ratio')
      }
    )
  
  rtpcr_export(d)
  rtpcr_export_nice(d)
  
  invisible(d)
  
}


read_time_series <- function(dataset, ranges = TS_RANGES) {
  
  ranges %>%
    purrr::map2(
      names(.),
      function(range, t) {
        read_and_preprocess(
          dataset = dataset,
          xlabels = XLABELS_1_2,
          range = range
        ) %>%
          dplyr::mutate(hours = readr::parse_number(t))
      }
    ) %>%
    dplyr::bind_rows()
  
}


read_and_preprocess <- function(
  dataset,
  xlabels = XLABELS_1_2,
  range = NULL
) {
  
  dataset %>%
    sprintf('AJ set %s.xlsx', .) %>%
    file.path() %>%
    {suppressMessages(readxl::read_excel(., sheet = 1L, range = range))} %>%
    rlang::set_names(
      dplyr::slice(., 1L) %>%
        as.character %>%
        make.unique(sep = '_') %>%
        stringr::str_replace_all(' ', '_')
    ) %>%
    dplyr::slice(2L:dplyr::n()) %>%
    dplyr::filter(!is.na(Targets)) %>%
    tidyr::pivot_longer(cols = c(-Targets, -Normalized)) %>%
    tidyr::separate(
      name,
      c('cp', 'variable', 'replicate'),
      sep = '_',
      fill = 'right'
    ) %>%
    dplyr::select(-cp) %>%
    dplyr::mutate(variable = stringr::str_to_upper(variable)) %>%
    tidyr::replace_na(list(replicate = "0")) %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    dplyr::rename(condition = Normalized, target = Targets) %>%
    tidyr::pivot_wider(
      names_from = variable,
      values_from = value
    ) %>%
    dplyr::filter(!is.na(TARG) & condition %in% xlabels) %>%
    dplyr::mutate(
      D_CT = TARG - REF,
      targ_by_ref = 2 ** -D_CT,
      control = condition %in% CONTROLS
    ) %>%
    dplyr::group_by(target, condition) %>%
    dplyr::mutate(
      SEM_D_CT = sd(D_CT) / sqrt(dplyr::n()),
      sem_targ_by_ref = sd(targ_by_ref) / sqrt(dplyr::n())
    ) %>%
    dplyr::ungroup() %>%
    {dplyr::inner_join(
      .,
      dplyr::filter(., control),
      by = c('target', 'replicate'),
      suffix = c('_d', '_c')
    )} %>%
    dplyr::filter(!control_d | condition_d == condition_c) %>%
    dplyr::filter(
      ifelse(
        condition_d %in% MET_CONTROL,
        condition_c == 'MET',
        condition_c == 'ET'
      )
    ) %>%
    dplyr::mutate(
      DD_CT = D_CT_d - D_CT_c,
      fold_induction = 2 ** -DD_CT
    ) %>%
    dplyr::group_by(target, condition_d) %>%
    dplyr::mutate(
      SEM_DD_CT = sd(DD_CT) / sqrt(dplyr::n()),
      sem_fold_induction = sd(fold_induction) / sqrt(dplyr::n()),
      pval_fold_induction =
        t.test(
          fold_induction,
          y = rep(1., 3L),
          paired = TRUE
        )$p.value,
      pval_targ_by_ref_d =
        t.test(
          targ_by_ref_d,
          y = targ_by_ref_c,
          paired = TRUE
        )$p.value,
      stars_fold_induction = stars_pval(pval_fold_induction),
      stars_targ_by_ref_d = stars_pval(pval_targ_by_ref_d),
      mean_fold_induction = mean(fold_induction),
      mean_targ_by_ref_d = mean(targ_by_ref_d)
    ) %>%
    dplyr::ungroup() %>%
    mutate(
      Receptor = ifelse(condition_d %in% MET_CONTROL, 'PXR', 'VDR'),
      condition_d = factor(condition_d, levels = xlabels, ordered = TRUE),
    ) %>%
    {`if`('hours' %in% names(.), ., mutate(., hours = 24))} %>%
    `attr<-`('idx', dataset)
  
}


rtpcr_barplot <- function(d, var, ylabel) {
  
  var <- rlang::enquo(var)
  var_chr <- rlang::quo_text(var)
  var_sem <- as.symbol(sprintf('sem_%s', var_chr))
  var_sta <- as.symbol(sprintf('stars_%s', var_chr))
  
  n_targets <-
    d %>%
    dplyr::pull(target) %>%
    dplyr::n_distinct()
  
  width <-
    n_targets %>%
    divide_by(2L) %>%
    ceiling %>%
    multiply_by(800L) %>%
    add(300L)
  
  height <- n_targets %>% is_less_than(3L) %>% `if`(1000L, 1800L)
  
  y_scale <-
    d %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(!!var := max(!!var)) %>%
    dplyr::summarize_all(first) %>%
    dplyr::pull(!!var) %>%
    {`if`(
      max(.) / min(.) > 100L,
      'free_y',
      'fixed'
    )}
  
  p <-
    ggplot2::ggplot(
      d %>%
        dplyr::group_by(target, condition_d) %>%
        dplyr::mutate(!!var := mean(!!var)) %>%
        dplyr::summarize_all(dplyr::first),
      ggplot2::aes(y = !!var, x = condition_d, fill = Receptor)
    ) +
    ggplot2::facet_wrap(
      ~target,
      nrow = 2L,
      scales = y_scale
    ) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = !!var - !!var_sem,
        ymax = !!var + !!var_sem
      ),
      width = .1
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = !!var_sta),
      color = 'red',
      size = 8,
      nudge_x = .2
    ) +
    ggplot2::theme_linedraw() +
    ggplot2::ggtitle('Expressions at transcript level under drug treatments') +
    ggplot2::xlab('Condition (drug)') +
    ggplot2::ylab(ylabel)
  
  dir.create(path='figures', showWarnings = FALSE)
  
  file.path(
    'figures',
    sprintf('%s_barplot_%s.png', var_chr, attr(d, 'idx'))
  ) %>%
    png(width = width, height = height, res = 300L)
  suppressWarnings(print(p))
  dev.off()
  
}


rtpcr_summary <- function(...) {
  
  d <-
    dplyr::bind_rows(...) %>%
    dplyr::group_by(target, condition_d) %>%
    dplyr::mutate(fold_induction = mean(fold_induction)) %>%
    dplyr::summarize_all(dplyr::first) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!condition_d %in% CONTROLS)
  
  p <-
    ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = condition_d,
        y = target,
        fill = log2(fold_induction),
        label = stars_fold_induction
      )
    ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = 'gray70', size = 10, nudge_y = -.2) +
    ggplot2::scale_fill_viridis_c(
      guide = guide_colorbar(title = 'Mean fold\ninduction\n(log2)')
    ) +
    ggplot2::ggtitle(
      'Fold induction of\nxenobiotic metabolism\ngenes by various compounds'
    ) +
    ggplot2::ylab('Xenobiotic metabolism genes') +
    ggplot2::xlab('Compounds') +
    ggplot2::theme_linedraw()
  
  dir.create(path='figures', showWarnings = FALSE)
  
  file.path('figures', 'fold_induction_summary.png') %>%
    png(width = 1000L, height = 1500L, res = 300L)
  suppressWarnings(print(p))
  dev.off()
  
}

rtpcr_summary_t <- function(...) {
  
  d <-
    list(...) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(target, condition_d, hours) %>%
    dplyr::mutate(fold_induction = mean(fold_induction)) %>%
    dplyr::summarize_all(dplyr::first) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(hours = as.character(hours)) %>%
    dplyr::filter(!condition_d %in% CONTROLS)
  
  p <-
    ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = hours,
        y = condition_d,
        fill = log2(fold_induction),
        label = stars_fold_induction
      )
    ) +
    ggplot2::facet_wrap(~target) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = 'gray70', size = 6, nudge_y = -.25) +
    ggplot2::scale_fill_viridis_c(
      guide = guide_colorbar(title = 'Mean fold\ninduction\n(log2)')
    ) +
    ggplot2::ggtitle(
      'Fold induction of xenobiotic metabolism\ngenes by various compounds'
    ) +
    ggplot2::ylab('Compounds') +
    ggplot2::xlab('Time (hours)') +
    ggplot2::theme_linedraw()
  
  dir.create(path='figures', showWarnings = FALSE)
  
  file.path('figures', 'fold_induction_summary_t.png') %>%
    png(width = 1200L, height = 1200L, res = 300L)
  suppressWarnings(print(p))
  dev.off()
  
}


rtpcr_times_plot <- function(d, var, ylabel) {
  
  var <- rlang::enquo(var)
  var_chr <- rlang::quo_text(var)
  var_sem <- as.symbol(sprintf('sem_%s', var_chr))
  var_sta <- as.symbol(sprintf('stars_%s', var_chr))
  
  n_targets <-
    d %>%
    dplyr::pull(target) %>%
    dplyr::n_distinct()
  
  width <-
    n_targets %>%
    divide_by(2L) %>%
    ceiling %>%
    multiply_by(800L) %>%
    add(300L)
  
  height <- n_targets %>% is_less_than(3L) %>% `if`(1000L, 1800L)
  
  y_scale <-
    d %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(!!var := max(!!var)) %>%
    dplyr::summarize_all(first) %>%
    dplyr::pull(!!var) %>%
    {`if`(
      max(.) / min(.) > 100L,
      'free_y',
      'fixed'
    )}
  
  d %<>%
    mutate(hours = as.character(hours)) %>%
    {`if`(var_chr == 'fold_induction', filter(., !control_d), .)}
  
  dodge <- ggplot2::position_dodge(.2)
  
  p <-
    ggplot2::ggplot(
      d %>%
        dplyr::group_by(target, condition_d, hours) %>%
        dplyr::mutate(!!var := mean(!!var)) %>%
        dplyr::summarize_all(dplyr::first),
      ggplot2::aes(
        y = !!var,
        x = hours,
        color = condition_d,
        group = condition_d
      )
    ) +
    ggplot2::facet_wrap(
      ~target,
      nrow = 2L,
      scales = y_scale
    ) +
    ggplot2::geom_line(
      position = dodge
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = !!var - !!var_sem,
        ymax = !!var + !!var_sem
      ),
      position = dodge,
      width = .4
    ) +
    ggplot2::geom_point(
      fill = 'white',
      shape = 21,
      position = dodge
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = !!var_sta),
      size = 8,
      nudge_x = .1,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_discrete(
      guide = ggplot2::guide_legend(title = 'Compound')
    ) +
    ggplot2::scale_x_discrete() +
    ggplot2::theme_linedraw() +
    ggplot2::ggtitle('Expressions at transcript level under drug treatments') +
    ggplot2::xlab('Time (hours)') +
    ggplot2::ylab(ylabel)
  
  dir.create(path='figures', showWarnings = FALSE)
  
  file.path(
    'figures',
    sprintf('%s_times_%s.png', var_chr, attr(d, 'idx'))
  ) %>%
    png(width = width, height = height, res = 300L)
  suppressWarnings(print(p))
  dev.off()
  
}


rtpcr_boxplot <- function(d, var, ylabel) {
  
  p <-
    ggplot2::ggplot(d, ggplot2::aes(y = fold_induction, x = condition_d)) +
    ggplot2::facet_wrap(~target, scales = 'free_y') +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(
      size = 3,
      color = 'red',
      width = .3,
      alpha = .5,
      shape = 18
    ) +
    ggplot2::theme_linedraw() +
    ggplot2::ggtitle('Expressions at transcript level under drug treatments') +
    ggplot2::xlab('Condition (drug)') +
    ggplot2::ylab('Fold induction')
  
  file.path('figures', 'fold_induction_boxplot.png') %>%
    png(width = 1600L, height = 1800L, res = 300L)
  print(p)
  dev.off()
  
}


stars_pval <- function(pval) {
  
  ifelse(pval > 0.0001, '***', '****') %>%
    ifelse(pval > 0.001, '**', .) %>%
    ifelse(pval > 0.01, '*', .) %>%
    ifelse(pval > 0.05, '', .)
  
}


rtpcr_export <- function(d, dataset = NULL) {
  
  dir.create('tables', showWarnings = FALSE)
  
  dataset %>%
    `%||%`(attr(d, 'idx')) %>%
    sprintf('rtpcr_set_%s.tsv', .) %>%
    file.path('tables', .) %>%
    readr::write_tsv(d, ., progress = FALSE)
  
}


rtpcr_export_nice <- function(d, dataset = NULL) {
  
  d %>%
    select(
      Gene = target,
      Condition = condition_d,
      Control = condition_c,
      Replicate = replicate,
      Time = hours,
      `Target/ref` = targ_by_ref_d,
      `Fold induction` = fold_induction,
      `P-value` = pval_fold_induction
    ) %>%
    `attr<-`('idx', sprintf('%s_summary', attr(., 'idx'))) %>%
    rtpcr_export(dataset = dataset)
  
}


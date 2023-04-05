
```
export_plots_2_png <- function(l, path = tempdir(), nrow = NULL, ncol = NULL, 
                               width = 12, height = 8) 
{
  len <- length(l)
  
  file_names <- here::here(path, paste0("fig_", 1:len, ".png"))
  file_names %>% 
    map2(l, function(x, y) {
      ggsave(filename = x, plot = y, width = width, height = height)
    })
  
  file_names
}

# combine png files -------------------------------------------------------

combine_png <- function(fig_paths, filename = "Plots_Combined%03d.png",
                        width = 15, height = 6, nrow = 2, ncol = 1,
                        title = NULL)
{
  plots <- lapply(fig_paths, function(x){
    img <- grDevices::as.raster(png::readPNG(x))
    grid::rasterGrob(img, interpolate = FALSE)
  })
  
  ggsave(filename = filename, width = width, height = height,
         gridExtra::marrangeGrob(
           grobs = plots, nrow = nrow, ncol = ncol, top = title
         ))
  
  fig_paths %>%
    purrr::map(fs::file_delete)
}




# plot flows by hpv -------------------------------------------------------



plot_flow_hpv <- function(data, title, count = TRUE, ymax) 
{
  
  df <- data %>% 
    filter(visit == "Screen") %>% 
    mutate(visit = fct_drop(visit)) 
  
  pvalues <- df %>% 
    rstatix::group_by(visit) %>% 
    rstatix::wilcox_test(x ~ hpv_result) %>% 
    rstatix::adjust_pvalue() %>% 
    rstatix::add_xy_position() %>%
    rstatix::add_significance()
  print(pvalues)
  
  # ds <- df %>% 
  #   group_by(visit) %>% 
  #   summarize(n = n()) %>% 
  #   ungroup() %>% 
  #   mutate(lbl = paste0(visit, " (", n, ")"))
  # lbl <- pull(ds, lbl)
  # names(lbl) <- pull(ds, visit)
  
  ds <- df %>%
    group_by(visit, hpv_result) %>%
    summarize(n = n()) %>%
    ungroup()
  
  p <- df %>% 
    ggplot(aes(x = hpv_result, y = x)) + 
    geom_boxplot(aes(fill = hpv_result), width = 0.4, alpha = .5) +
    scale_fill_manual(values = c("red", "blue")) +
    geom_jitter(aes(fill = hpv_result), position = position_jitterdodge(0.4), 
                size=2, shape=21, stroke = 1) +
    scale_color_manual(values = c("red", "blue")) + 
    theme_solid_frame 
  
  if (count) {
    ylab <- "Log10"
    p <- p + 
      geom_text(data = ds, aes(
        y = 0.5, x = hpv_result, label = paste0("n=", n), color = hpv_result
      )) + 
      scale_y_log10(
        limits = c(-1, ymax),
        breaks = c(0, 1, 100, 1e4, 1e6),
        label = scales::comma
      )
  } else {
    ylab <- "Percentage"
    p <- p + 
      geom_text(data = ds, aes(
        y = -10, x = hpv_result, label = paste0("n=", n), color = hpv_result
      )) + 
      scale_y_continuous(
        limits = c(-10, ymax),
        breaks = seq(0, 100, 25)
      )
  }
  p + 
    labs(
      title = title, y = ylab, x = "HPV Status", color = "HPV", fill = "HPV"
    ) +
    theme(legend.position = "none")
}

plot_flow_visit <- function(data, title, count = TRUE, ymax) 
{
  df <- data %>% 
    filter(hpv_result == "Positive") %>% 
    mutate(hpv_result = fct_drop(hpv_result)) 
  
  pvalues <- df %>% 
    rstatix::group_by(hpv_result) %>% 
    rstatix::wilcox_test(x ~ visit, ref.group = "Enrol") %>% 
    filter(p.adj < 0.05) %>%
    rstatix::add_xy_position()
  print(pvalues)
  
  ds <- data %>% 
    group_by(visit) %>% 
    summarize(n = n()) %>% 
    ungroup() %>% 
    mutate(lbl = paste0(visit, " (", n, ")"))
  
  p <- df %>% 
    ggplot(aes(x = visit, y = x)) +
    geom_boxplot(aes(fill = visit), width = 0.4, alpha = .5) +
    geom_jitter(aes(fill = visit), width = 0.1, 
                size=1.5, shape=21, stroke = 1) +
    theme_solid_frame
  
  if (count) {
    ylab <- "Log10"
    p <- p +
      geom_text(data = ds, aes(
        y = 0.5, x = visit, label = paste0("n=", n), color = visit
      )) + 
      scale_y_log10(
        limits = c(-1, ymax),
        breaks = c(0, 1, 100, 1e4, 1e6),
        label = scales::comma
      )
    if (nrow(pvalues) > 0) p <- p + 
      ggprism::add_pvalue(
        pvalues, step.increase = 0.1
      )
  } else {
    ylab <- "Percentage"
    p <- p + 
      geom_text(data = ds, aes(
        y = -10, x = visit, label = paste0("n=", n), color = visit
      )) + 
      scale_y_continuous(
        limits = c(-10, ymax),
        breaks = seq(0, 100, 25)
      )
    if (nrow(pvalues) > 0) p <- p +
      ggprism::add_pvalue(pvalues)
  }
  p + 
    labs(
      title = NULL, y = NULL, x = "Visits (HPV+)", color = "HPV", fill = "HPV"
    ) +
    theme(legend.position = "none")
}





plot_flow_ct <- function(data, x, title)
{
  data %>%
    mutate(x = {{ x }}) %>% 
    select(pid, visit, x, hr_hpv_16:hr_p5) %>% 
    pivot_longer(hr_hpv_16:hr_p5, names_to = "ct") %>% 
    filter(!is.na(value)) %>% 
    ggplot(aes(x = x, y = value, color = visit)) + 
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_y_log10(label = scales::comma) +
    scale_x_log10(label = scales::math_format(format = log10)) +
    facet_grid(visit ~ ct, labeller = label_both) + 
    labs(
      title = title, 
      x = NULL, y = "ct values", color = "Visit"
    ) + 
    theme_solid_frame +
    theme(legend.position = "bottom")
}

 ``` 

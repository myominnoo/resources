

theme_scale_y_log10 <- scale_y_log10(
  limits = c(0, 1e8), 
  breaks = c(1, 100, 1e4, 1e6),
  label = scales::comma
)

theme_scale_y_log10_fixed <- scale_y_log10(
  breaks = c(1, 100, 1e4, 1e6),
  label = scales::math_format(format = log10)
)

theme_solid_frame <- theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill="white")
  )
theme_rotate_x_axis_label <- theme(
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
)

theme_rotate_x_axis_label_45 <- theme(
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)


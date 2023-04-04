

# themes ------------------------------------------------------------------

theme_scale_y_log10 <- scale_y_log10(breaks = c(1, 100, 1e4),
																		 label = scales::comma)

theme_scale_y_log10_fixed <- scale_y_log10(
	limits = c(1, 1e4),
	breaks = c(1, 100, 1e4),
	label = scales::comma
)

theme_solid_frame <- theme_bw() +
	theme(
		text = element_text(size = 12),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background=element_rect(fill="white")
	)

theme_scale_x_fixed <- scale_x_continuous(
	limits = c(-100, 650), breaks = seq(-100, 650, 100)
)

theme_rotate_x_axis_label <- theme(
	axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
)


# Figure functions --------------------------------------------------------

plot_fig1.1 <- function(data, y, title = NULL,
												xaxis = "Time since initial dose (days)",
												legend.position = "none") {
	data %>%
		ggplot(aes(x = time_dose_1, y = {{ y }}, group = pid,
							 color = group, shape = group)) +
		geom_point() +
		geom_line(alpha = 0.15) +
		geom_hline(yintercept = 1, alpha = 0.5, linetype = "dotted") +
		geom_vline(xintercept = 0, alpha = 0.5, linetype = "dotted") +
		theme_scale_y_log10_fixed +
		labs(
			title = title,
			x = xaxis, y = "BAU/ml",
			color = "Cohort", shape = "Cohort"
		) +
		theme_solid_frame +
		theme(legend.position = legend.position)
}

plot_fig1.2 <- function(data, y, title = NULL,
												xaxis = "Time since initial dose (days)",
												legend.position = "none") {
	data %>%
		ggplot(aes(x = time_dose_1, y = {{ y }}, group = pid,
							 color = pre_vax_infection, shape = pre_vax_infection)) +
		geom_point() +
		geom_line(alpha = 0.15) +
		geom_hline(yintercept = 1, alpha = 0.5, linetype = "dotted") +
		geom_vline(xintercept = 0, alpha = 0.5, linetype = "dotted") +
		theme_scale_y_log10_fixed +
		labs(
			title = title,
			x = xaxis, y = "BAU/ml",
			color = "Cohort", shape = "Cohort"
		) +
		theme_solid_frame +
		theme(legend.position = legend.position)
}

plot_fig1.3 <- function(data, x, title) {
	df <- data %>%
		mutate(serology_rbd = {{ x }})
	pvalues <- df %>%
		mutate(serology_rbd = log10(serology_rbd)) %>%
		# group_by(time_since_dose_x_cat, pre_vax_infection) %>%
		# summarize(n = n()) %>%
		# filter(n > 5)
		# View
		filter(time_since_dose_x_cat %in% c(
			"D1-D2: ≤14 days", "D2-D3: 2-3 mnth"
		)) %>%
		rstatix::group_by(time_since_dose_x_cat) %>%
		rstatix::wilcox_test(serology_rbd ~ pre_vax_infection) %>%
		rstatix::add_xy_position() %>%
		rstatix::add_significance()

	## figure 2
	df %>%
		filter(!time_since_dose_x_cat %in% c(
			"Post-D5: ≤1 mnth", "Post-D5: 1+ mnth"
		)) %>%
		ggplot(aes(x = pre_vax_infection, y = serology_rbd,
							 color = pre_vax_infection)) +
		geom_boxplot(color = "black") +
		geom_jitter(size = 1,
								# position = position_jitterdodge(0.2),
								alpha = 0.5) +
		scale_color_manual(values = c("blue", "red")) +
		theme_scale_y_log10 +
		labs(
			x = "Prior COVID-19 infection", y = "BAU/ml"
		) +
		facet_wrap(~ time_since_dose_x_cat, nrow = 1) +
		theme_solid_frame +
		theme(legend.position = "none",
					strip.background=element_rect(fill="white")) +
		ggpubr::stat_pvalue_manual(pvalues)
}

plot_fig1.4 <- function(data, x) {
	df <- data %>%
		mutate(serology_rbd = {{ x }})
	pvalues <- df %>%
		filter(time_point3 %in% c("Pre-D1", "D1-D2", "D2-D3", "D3-D4")) %>%
		mutate(serology_rbd = log10(serology_rbd)) %>%
		rstatix::group_by(time_since_dose_x_cat) %>%
		rstatix::tukey_hsd(serology_rbd ~ group) %>%
		rstatix::add_xy_position() %>%
		rstatix::add_significance() %>%
		filter(p.adj < 0.05)
	## figure 2
	df %>%
		filter(!time_since_dose_x_cat %in% c(
			"Post-D5: ≤1 mnth", "Post-D5: 1+ mnth"
		)) %>%
		ggplot(aes(x = group, y = serology_rbd,
							 color = group)) +
		geom_boxplot(color = "black") +
		geom_jitter(size = 1,
								# position = position_jitterdodge(0.2),
								alpha = 0.5) +
		theme_scale_y_log10 +
		labs(y = "BAU/ml") +
		facet_wrap(~ time_since_dose_x_cat, nrow = 1) +
		theme_solid_frame +
		theme_rotate_x_axis_label +
		theme(legend.position = "none",
					strip.background=element_rect(fill="white")) +
		ggpubr::stat_pvalue_manual(pvalues)

}





plot_fig2.1 <- function(data, x, y, title = NULL,
											 legend.position = "none",
											 color = "Cohort")
{
	data %>%
		ggplot(aes(x = {{ x }}, y = {{ y }},
							 color = {{ x }})) +
		geom_boxplot(color = "black") +
		geom_jitter(size = 1) +
		theme_scale_y_log10 +
		labs(
			title = title,
			x = NULL, y = "BAU/ml",
			color = color, shape = color
		) +
		theme_solid_frame +
		theme(legend.position = legend.position,
					strip.background=element_rect(fill="white")) +
		theme_rotate_x_axis_label
}


combine_plot_fig2.1 <- function(data, x, title = NULL)
{
	data <- data %>%
		bind_rows(
			ab_mod %>%
				mutate(group = "All")
		) %>%
		mutate(serology_rbd = {{ x }},
					 group = fct_expand(group, "All"))

	pval1 <- data %>%
		filter(time_since_dose_x_cat %in% c(
			"Pre-D1", "D1-D2: ≤14 days", "D1-D2: 14+ days"
		)) %>%
		mutate(serology_rbd = log10(serology_rbd)) %>%
		rstatix::group_by(group) %>%
		rstatix::wilcox_test(serology_rbd ~ time_since_dose_x_cat,
												 ref.group = "D1-D2: ≤14 days") %>%
		rstatix::add_xy_position()

	pval2 <- data %>%
		filter(group != "GER") %>%
		filter(time_since_dose_x_cat %in% c(
			"D1-D2: 14+ days",
			"D2-D3: ≤1 mnth", "D2-D3: 2-3 mnth", "D2-D3: 4-6 mnth", "D2-D3: 6+ mnth"
		)) %>%
		mutate(serology_rbd = log10(serology_rbd)) %>%
		rstatix::group_by(group) %>%
		rstatix::wilcox_test(serology_rbd ~ time_since_dose_x_cat,
												 ref.group = "D2-D3: ≤1 mnth") %>%
		rstatix::add_xy_position()

	pval3 <- data %>%
		filter(group != "CKD") %>%
		filter(time_since_dose_x_cat %in% c(
			"D2-D3: 6+ mnth",
			"D3-D4: ≤1 mnth", "D3-D4: 2-3 mnth", "D3-D4: 4-6 mnth", "D3-D4: 6+ mnth"
		)) %>%
		mutate(serology_rbd = log10(serology_rbd)) %>%
		rstatix::group_by(group) %>%
		rstatix::wilcox_test(serology_rbd ~ time_since_dose_x_cat,
												 ref.group = "D3-D4: ≤1 mnth") %>%
		rstatix::add_xy_position()

	pval3_ckd <- data %>%
		filter(group == "CKD") %>%
		filter(time_since_dose_x_cat %in% c(
			"D2-D3: 6+ mnth",
			"D3-D4: ≤1 mnth", "D3-D4: 2-3 mnth", "D3-D4: 4-6 mnth", "D3-D4: 6+ mnth"
		)) %>%
		mutate(serology_rbd = log10(serology_rbd)) %>%
		rstatix::group_by(group) %>%
		rstatix::wilcox_test(serology_rbd ~ time_since_dose_x_cat,
												 ref.group = "D3-D4: 4-6 mnth") %>%
		rstatix::add_xy_position()

	pvalues <- bind_rows(pval1, pval2, pval3, pval3_ckd) %>%
		filter(p.adj < 0.05)
	pvalues

	data %>%
		plot_fig2.1(time_since_dose_x_cat, serology_rbd, title = title) +
		facet_wrap(~ group, nrow = 1) +
		ggpubr::stat_pvalue_manual(pvalues)
}



plot_fig2.2 <- function(data, x, filter_time = "D2-D3", title = "(A) D2-D3",
												yline = 4.2, ylim = 5e4, ytext = 3e4)
{
	lmm_ds <- data %>%
		mutate(
			serology_rbd = {{ x }},
			time_since_dose_x = as.numeric(time_since_dose_x),
			time_since_dose_x_mnth = round(time_since_dose_x / 30.6, 0)
		) %>%
		filter(time_point3 %in% filter_time) %>%
		mutate(
			time_since_dose_x_mnth = as.numeric(time_since_dose_x_mnth),
			time_group = case_when(
				time_since_dose_x <= 30 ~ 30,
				time_since_dose_x > 30 & time_since_dose_x <= 90 ~ 90,
				time_since_dose_x > 90 & time_since_dose_x <= 180 ~ 180,
				time_since_dose_x > 180 ~ 270
			),
			time_since_dose_x_mnth2 = case_when(
				time_since_dose_x_mnth <= 1 ~ 1,
				time_since_dose_x_mnth > 1 & time_since_dose_x_mnth <= 3 ~ 3,
				time_since_dose_x_mnth > 3 & time_since_dose_x_mnth <= 6 ~ 6,
				time_since_dose_x_mnth > 6 ~ max(.$time_since_dose_x_mnth)
			)
		)

	if (filter_time == "D2-D3") {
		lmm_ds <- lmm_ds %>%
			mutate(
				time_since_dose_x_mnth_cutoff = ifelse(
					time_since_dose_x <= median(.$time_since_dose_x), 0, 1
				)
			)
	} else {
		lmm_ds <- lmm_ds %>%
			mutate(
				time_since_dose_x_mnth_cutoff = ifelse(
					time_since_dose_x <= 90, 0, 1
				)
			)
	}

	## months cut-off
	p1 <- lmm_ds %>%
		ggplot(aes(x = time_since_dose_x, y = serology_rbd, group = pid,
							 color = group, shape = group)) +
		geom_point() +
		geom_line(alpha = 0.15) +
		geom_hline(yintercept = 1, linetype = "dotted") +
		geom_vline(xintercept = 0, linetype = "dotted") +
		geom_vline(xintercept = 30*3, linetype = "dotted", color = "red") +
		geom_vline(xintercept = 30*6, linetype = "dotted", color = "red") +
		geom_vline(xintercept = 30.5*12, linetype = "dotted", color = "red") +
		labs(x = "Days", y = "BAU/ml") +
		scale_color_manual(values = scales::hue_pal()(3)) +
		scale_y_log10(
			limits = c(1, ylim),
			breaks = c(1, 100, 1e4),
			label = scales::comma
		) +
		scale_x_continuous(limits = c(0, 540), breaks = seq(0, 550, 90)) +
		theme_solid_frame +
		theme(legend.position = "none")
	p1

	ph1 <- lmerTest::lmer(serology_rbd ~ time_since_dose_x_mnth +
													(1|pid), data = lmm_ds)
	summary(ph1)
	cc <- coef(summary(ph1))
	decay_all <- floor(abs(cc[grepl("time_since_dose_x", rownames(cc)),][1] / 12))
	decay_all_pval <- cc[grepl("time_since_dose_x", rownames(cc)),][5] %>%
		data.frame(p = .) %>%
		rstatix::add_significance() %>%
		pull(p.signif)


	ph2 <- lmerTest::lmer(serology_rbd ~ time_since_dose_x_mnth +
													time_since_dose_x_mnth_cutoff +
													(1|pid), data = lmm_ds)
	cc <- coef(summary(ph2))
	decay_t1 <- floor(abs(cc[grepl("time_since_dose_x_mnth",
																 rownames(cc)),][1] / 12))
	decay_t2 <- floor(abs(cc[grepl("time_since_dose_x_mnth_cutoff",
																 rownames(cc)),][1] / 12))
	decay_diff_pval <- cc[grepl("time_since_dose_x_mnth_cutoff",
															rownames(cc)),][5] %>%
		data.frame(p = .) %>%
		rstatix::add_significance() %>%
		pull(p.signif)

	p2 <- bind_rows(
		lmm_ds %>%
			group_by(time_since_dose_x_mnth2) %>%
			summarise(mean = mean(serology_rbd, na.rm = TRUE),
								sd = sd(serology_rbd, na.rm = TRUE)) %>%
			ungroup() %>%
			mutate(group = "All"),
		lmm_ds %>%
			group_by(group, time_since_dose_x_mnth2) %>%
			summarise(mean = mean(serology_rbd, na.rm = TRUE),
								sd = sd(serology_rbd, na.rm = TRUE)) %>%
			ungroup()
	) %>%
		mutate(
			group = factor(group, levels = c("CKD", "NML1", "NML2", "All")),
			n = n(),
			err = qt(0.975, df = n-1)*sd/sqrt(n),
			lower = mean - err,
			upper = mean + err,
			across(lower:upper, ~ ifelse(is.na(.x) | .x < 0, 10, .x))
		) %>%
		ggplot(aes(x = time_since_dose_x_mnth2, y = mean,
							 ymin = lower, ymax = upper, color = group, shape = group)) +
		geom_point(position = position_dodge(width = 1), size = 2) +
		geom_errorbar(position = position_dodge(width = 1), width = 1) +
		geom_line(aes(linetype = group), position = position_dodge(width = 1)) +
		geom_hline(yintercept = 1, linetype = "dotted") +
		geom_vline(xintercept = 0, linetype = "dotted") +
		geom_segment(x = 0, y = yline, xend = 2.9, yend = yline,
								 color = "black") +
		geom_segment(x = 3.1, y = yline, xend = max(lmm_ds$time_since_dose_x_mnth),
								 yend = yline,
								 color = "black") +
		annotate(geom="text", x = 1.5, y = ytext, label= paste0(
			"t½ = ", decay_t1, " days"
		), color="black") +
		annotate(geom="text", x = 9, y = ytext, label= paste0(
			"t½ = ", decay_t2, " days"
		), color="black") +
		annotate(geom="text", x = 15, y = 5, label= paste0(
			"t½ = ", decay_all, " days, (p = ", decay_all_pval, ")"
		), color="black") +
		annotate(geom="text", x = 15, y = 3, label= paste0(
			"Δ decay rate = ", decay_diff_pval
		), color="black") +
		labs(
			x = "Months", y = "BAU/ml",
			color = "Cohort", shape = "Cohort", linetype = "Cohort"
		) +
		scale_color_manual(values = c(scales::hue_pal()(3), "#C77CFF")) +
		scale_y_log10(
			limits = c(1, ylim),
			breaks = c(1, 100, 1e4),
			label = scales::comma
		) +
		scale_x_continuous(limits = c(0, 18), breaks = seq(0, 20, 3)) +
		theme_solid_frame +
		theme(legend.position = "right")

	p1 + p2 +
		plot_layout(design = "AAABBB") +
		plot_annotation(title = title)
}




plot_fig3.1 <- function(data, title)
{
	df <- data %>%
		filter(!(time_point3 %in% c("D4-D5", "Post-D5"))) %>%
		pivot_longer(-c(1, 2, 3)) %>%
		separate(name, into = c("prnt", "type")) %>%
		mutate(type = toupper(type),
					 type = factor(type),
					 type = relevel(type, ref = "WT"))

	pvalues <- df %>%
		filter(time_point3 %in% c("D1-D2", "D2-D3")) %>%
		# mutate(value = log10(value)) %>%
		rstatix::group_by(time_point3) %>%
		rstatix::wilcox_test(value ~ type) %>%
		filter(p.adj < 0.05) %>%
		rstatix::add_xy_position() %>%
		mutate(y.position = seq(3.0, 5, 0.2)[1:n()])

	df %>%
		ggplot(aes(x = type, y = value, color = type)) +
		geom_boxplot(color = "black") +
		geom_jitter() +
		facet_wrap(~ time_point3, nrow = 1, ncol = 6) +
		labs(title = title, x = NULL, y = "BAU/ml") +
		theme_scale_y_log10 +
		theme_solid_frame +
		theme(legend.position = "none") +
		ggpubr::stat_pvalue_manual(pvalues)
}



plot_fig3.2 <- function(data, title)
{
	df <- data %>%
		select(pid, time_since_dose_x, time_point3,
					 lentivirus_neutralization_ic50_wt:lentivirus_neutralization_ic50_omicron_ba_4_5) %>%
		setNames(
			names(.) %>%
				gsub("lentivirus_neutralization_ic50_", "", .) %>%
				gsub("ba_1", "BA1", .) %>%
				gsub("ba_2", "BA2", .) %>%
				gsub("ba_4_5", "BA45", .)
		) %>%
		filter(!(time_point3 %in% c("D4-D5", "Post-D5"))) %>%
		pivot_longer(-c(1, 2, 3)) %>%
		separate(name, into = c("type", "voc")) %>%
		mutate(type = toupper(type),
					 type = ifelse(type == "OMICRON", voc, type),
					 type = factor(type),
					 type = relevel(type, ref = "WT"))

	pvalues1 <- df %>%
		filter(time_point3 %in% c("D2-D3")) %>%
		# mutate(value = log10(value)) %>%
		rstatix::group_by(time_point3) %>%
		rstatix::wilcox_test(value ~ type) %>%
		filter(p.adj < 0.05) %>%
		rstatix::add_xy_position() %>%
		mutate(y.position = seq(5, 6.5, 0.2)[1:n()])
	pvalues2 <- df %>%
		filter(time_point3 %in% c("D3-D4")) %>%
		# mutate(value = log10(value)) %>%
		rstatix::group_by(time_point3) %>%
		rstatix::wilcox_test(value ~ type) %>%
		filter(p.adj < 0.05) %>%
		rstatix::add_xy_position() %>%
		mutate(y.position = seq(5, 6.5, 0.2)[1:n()])
	pvalues <- bind_rows(pvalues1, pvalues2)

	df %>%
		ggplot(aes(x = type, y = value, color = type)) +
		geom_boxplot(color = "black") +
		geom_jitter() +
		facet_wrap(~ time_point3, nrow = 1, ncol = 6) +
		labs(title = title, x = NULL, y = "BAU/ml") +
		theme_scale_y_log10 +
		theme_solid_frame +
		theme(legend.position = "none") +
		ggpubr::stat_pvalue_manual(pvalues)
}



plot_fig3.3 <- function(data, title)
{
	df <- data %>%
		select(pid, time_since_dose_x, time_point3,
					 msd_neutralization_ic50_wt:msd_neutralization_ic50_omicron_ba_4_5) %>%
		setNames(
			names(.) %>%
				gsub("msd_neutralization_ic50_", "", .) %>%
				gsub("ba_1", "BA1", .) %>%
				gsub("ba_2", "BA2", .) %>%
				gsub("ba_4_5", "BA45", .)
		) %>%
		filter(!(time_point3 %in% c("D4-D5", "Post-D5"))) %>%
		pivot_longer(-c(1, 2, 3)) %>%
		separate(name, into = c("type", "voc")) %>%
		mutate(type = toupper(type),
					 type = ifelse(type == "OMICRON", voc, type),
					 type = factor(type),
					 type = relevel(type, ref = "WT"))
	df %>%
		ggplot(aes(x = type, y = value, color = type)) +
		geom_boxplot(color = "black") +
		geom_jitter() +
		facet_wrap(~ time_point3, nrow = 1, ncol = 6) +
		labs(title = title, x = NULL, y = "BAU/ml") +
		theme_scale_y_log10 +
		theme_solid_frame +
		theme(legend.position = "none")
}


plot_fig4.1 <- function(data, data2, x, title = NULL)
{

	df <- data %>%
		arrange(pid, desc(date_collected)) %>%
		distinct(pid, .keep_all = TRUE) %>%
		mutate(bt = "Pre") %>%
		bind_rows(data2) %>%
		# select(pid, group, date_collected, time_since_dose_x, time_point3, bt,
		# 			 serology_rbd) %>%
		mutate(serology_rbd = {{ x }}) %>%
		mutate(bt = ifelse(is.na(bt), "Post", bt),
					 bt = factor(bt, levels = c("Pre", "Post"))) %>%
		filter(!(time_point3 %in% c("D1-D2", "Post-D5")))

	pvalues <- df %>%
		rstatix::group_by(time_point3) %>%
		rstatix::wilcox_test(serology_rbd ~ bt, p.adjust.method = "BH") %>%
		rstatix::add_significance() %>%
		filter(p < 0.05) %>%
		rstatix::add_xy_position() %>%
		mutate(y.position = log10(y.position))
#
# 	df %>%
# 		ggplot(aes(x = bt, y = serology_rbd, color = bt)) +
# 		geom_boxplot(color = "black") +
# 		geom_jitter() +
# 		facet_wrap(~ time_point3, nrow = 1, ncol = 6) +
# 		labs(title = title, x = NULL, y = "BAU/ml") +
# 		scale_y_log10(
# 			limits = c(1, 1e5),
# 			breaks = c(1, 100, 1e4),
# 			label = scales::comma
# 		) +
# 		theme_solid_frame +
# 		theme(legend.position = "none") +
# 		ggpubr::stat_pvalue_manual(pvalues)
	df %>%
		ggplot(aes(x = bt, y = serology_rbd)) +
		geom_boxplot(aes(fill = bt), width = 0.7, alpha = .5) +
		geom_jitter(aes(fill = bt), position = position_jitterdodge(0.9),
								size=2, shape=21, stroke = 1) +
		scale_fill_manual(values = c("blue", "red")) +
		facet_wrap(~ time_point3, nrow = 1, ncol = 6) +
		labs(title = title, x = "Breakthrough infection", y = "BAU/ml") +
		scale_y_log10(
			limits = c(1, 1e5),
			breaks = c(1, 100, 1e4),
			label = scales::comma
		) +
		theme_solid_frame +
		theme(legend.position = "none") +
		ggpubr::stat_pvalue_manual(pvalues)
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
}







# calculate pvalues -------------------------------------------------------


calc_pvalues_tukey <- function(data, time_filter, xvar, yvar, ref.group,
															 step.increase = 0.12)
{
	data %>%
		filter(time_since_dose_x_cat %in% time_filter) %>%
		mutate(xvar = log10({{ xvar }}),
					 yvar = {{ yvar }}) %>%
		rstatix::tukey_hsd(xvar ~ yvar,
											 ref.group = ref.group) %>%
		rstatix::add_xy_position(step.increase = step.increase)
}


calc_pvalues <- function(data, time_filter, xvar, yvar, zvar, ref.group,
												 step.increase = 0.12)
{
	data %>%
		filter(time_since_dose_x_cat %in% time_filter) %>%
		mutate(x = log10({{ xvar }}),
					 y = {{ yvar }},
					 z = {{ zvar }}) %>%
		rstatix::group_by(z) %>%
		rstatix::wilcox_test(x ~ y, ref.group = ref.group) %>%
		rstatix::add_xy_position(step.increase = step.increase)
}


  

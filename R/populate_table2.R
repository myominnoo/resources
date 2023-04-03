populate_table1 <- function(data, by,
														byname = "Example",
														missing = c("no", "ifany", "always"),
														percent = c("column", "row"),
														add_N = TRUE,
														add_total = TRUE,
														add_pvalue = TRUE) {
	tbl <- data %>%
		gtsummary::tbl_summary(
			by = {{by}},
			statistic = list(gtsummary::all_continuous() ~ "{median} ({p25} - {p75})",
											 gtsummary::all_categorical() ~ "{n} ({p})"),
			digits = list(
				gtsummary::all_continuous() ~ 1,
				gtsummary::all_categorical() ~ c(0, 1)
			),
			missing = missing[1],
			missing_text = "Missing",
			percent = percent[1]
		)

	header <- dplyr::pull(gtsummary::show_header_names(tbl), "column")
	header <- header[grepl("stat_", header) & header != "stat_0"]
	tbl <- tbl %>%
		gtsummary::modify_spanning_header(header ~ sprintf("**%s**", byname))

	if (add_N) tbl <- tbl %>% gtsummary::add_n()
	if (add_total) {
		tbl <- tbl %>%
			gtsummary::add_overall(
				last = FALSE,
				statistic = list(gtsummary::all_continuous() ~ "{median} ({p25} - {p75})",
												 gtsummary::all_categorical() ~ "{n} ({p})"),
				digits = list(
					gtsummary::all_continuous() ~ 1,
					gtsummary::all_categorical() ~ c(0, 1)
				)
			)
	}
	if (add_pvalue) {
		tbl <- tbl %>%
			# rounding p-values to 3 decimal places
			gtsummary::add_p(
				pvalue_fun = function(x) gtsummary::style_number(x, digits = 3)
			)
	}
	tbl
}

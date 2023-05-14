
## inspired by someone from linkedIn.

library(tidyverse)
library(extrafont)

palmerpenguins::penguins |> 
	drop_na() |> 
	ggplot(aes(x = bill_length_mm, y = flipper_length_mm)) + 
	ggforce::geom_mark_ellipse(
		aes(fill = species, 
				fill = after_scale(colorspace::lighten(fill, .5)), 
				color = after_scale(colorspace::darken(fill, .5)), 
				label = species), 
		label.family = "Comic Sans MS", 
		label.fontsize = 12, 
		alpha = .2
	) + 
	geom_point(
		aes(fill = species, 
				fill = after_scale(colorspace::lighten(fill, .5)), 
				color = after_scale(colorspace::darken(fill, .5))), 
		shape = 21, 
		size = 2.2, 
		stroke = .75
	) + 
	scale_x_continuous(
		expand = c(.01, .01), 
		limits = c(30, 65), 
	) + 
	scale_y_continuous(
		expand = c(.01, .01), 
		limits = c(160, 260)
	) +
	labs(
		x = "Bill length (mm)", 
		y = "Flipper length (mm)", 
		title = "Bill versus flipper dimensions of brush-tailed penguins", 
		subtitle = "Scatterplot of bill versus flipper depth (colored ellipse marks show species grouping)", 
		caption = "Source: palmerpenguins-package {palmerpenguins}"
	) + 
	theme_light(base_size = 12, base_family = "Comic Sans MS") + 
	theme(
		legend.position = "none", 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		plot.title.position = "plot", 
		plot.title = element_text(face = "bold", size = 18), 
		plot.subtitle = element_text(face = "bold", color = "grey40", 
																 margin = margin(b = 15)), 
		plot.caption.position = "plot", 
		plot.caption = element_text(face = "bold.italic", color = "grey40"), 
		axis.title = element_text(face = "bold", color = "grey40"), 
		plot.margin = margin(rep(20, 4))
	)
ggsave("palmerpenguins-ecllipse.png", 
			 width = 8, height = 7, dpi = 300)	



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#                 HPV Microbiome Data Analysis
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Author: Myo Minn Oo
# Start Date: 30-06-2022


# prepare workspace -------------------------------------------------------

rm(list = ls())


## load packages 
library(tidyverse)
library(magrittr)
library(tidyverse)
library(janitor)
library(gtsummary)
library(flextable)
library(sjlabelled)

# theme to always run add_p() and bold_labels() after tbl_summary()
list(
  "tbl_summary-fn:addnl-fn-to-run" = 
    function(x) {
      bold_labels(x)  # bold labels and return table
    }
) %>%
  set_gtsummary_theme()

# import data
data.raw <- readxl::read_excel("data/MSM_HPV_cleaner_22june2022.xlsx")

# set high risk and low risk hpv types
hr_hpv_types <- c(31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 
                  26, 53, 73, 82, 66, 68)
lr_hpv_types <- c(6, 11, 30, 67, 70)

# clean and process data --------------------------------------------------

# what does cd4pc mean?
# HIV_code1, T.ARV, ARV, HIV_code2: same values different columns
# what does NA or missing means here? unknown?
# can we assume ART- to be ART naive and HIV+ to be on ART?


## subset variables for hpv, hiv and socio-demographics 
data.sd <- data.raw %>%
  # remove immunological and microbiological markers
  select(PID:HIV_code2, STI_TEST:O.INFO, HPV6:Number_of_HPV_types) %>% 
  select(PID, NMLTypingResult, HPV_pos, HIV, 
         B.YEAR,	AGE, AGE.F.SEX,	M.STATUS,	EDU.LEVEL, 
         CD4pc, CD4C, VL, HIV_code1, PREP, PEP,
         STI_TEST, STI_RESULT2,	STI.RESULT, 
         MSM.STATUS.x, S.WORKER, Y.S.PROST, D.PROST,
         REGULAR_partner, CASUAL_partner, P.P, ORAL_sex_reg, 
         ORAL_sex_cas, ORAL_sex_pay, INSERT_sex_reg, INSERT_sex_cas,
         INSERT_sex_pay, RECEPT_sex_reg, RECEPT_sex_cas, RECEPT_sex_pay,
         CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING, DRUG.HIST,
         C.M.PROB, O.INFO, 
         # HPV6:HPV97,
         Number_of_HPV_types
         ) 

## get NML typing result 
nml <- data.sd %>%
  select(PID, NMLTypingResult) %>% 
  filter(!grepl("nega|incon", NMLTypingResult)) %>%
  {
    nmax <- max(stringr::str_count(.$NMLTypingResult, ","), na.rm = TRUE) + 1
    separate(., NMLTypingResult, paste0("hpv.", seq_len(nmax)), 
             sep = ",", fill = "right")
  } %>% 
  pivot_longer(-PID) %>%
  select(-name) %>% 
  filter(!is.na(value)) %>% 
  mutate(value = as.numeric(value), 
         name = paste0("hpv_", value)) %>% 
  arrange(value) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(num_hpv = rowSums(!across(contains("hpv"), is.na)), 
         hr_hpv = rowSums(!across(num_range("hpv_", hr_hpv_types), is.na)),
         lr_hpv = rowSums(!across(num_range("hpv_", lr_hpv_types), is.na))) %>%
  unite(num_hpv_all, contains("hpv_"), sep = "+", na.rm = TRUE, remove = FALSE) %>% 
  unite(hr_hpv_all, num_range("hpv_", hr_hpv_types), sep = "+",
        na.rm = TRUE, remove = FALSE) %>% 
  unite(lr_hpv_all, num_range("hpv_", lr_hpv_types), sep = "+",
        na.rm = TRUE, remove = FALSE) %>% 
  select(PID, num_hpv_all, num_hpv, hr_hpv, hr_hpv_all, 
         lr_hpv_all, lr_hpv, hpv_16, hpv_18) %>%
  mutate(hpv_16 = ifelse(is.na(hpv_16), 0, 1), 
         hpv_18 = ifelse(is.na(hpv_18), 0, 1), 
         hpv_all = case_when(
           hpv_16 > 0 | hpv_18 > 0 ~ "HPV 16+18", 
           hr_hpv > 0 ~ "HR HPV", 
           lr_hpv > 0 ~ "LR HPV", 
           TRUE ~ "Other HPV"
         )) 

# combine data.sd and nml 
hpv <- data.sd %>% 
  mutate(
    # assume NA in HPV_pos as negative 
    HPV_pos = ifelse(is.na(HPV_pos), 0, HPV_pos),
    HPV_pos_cat = case_when(
      HPV_pos == 0 ~ "No", 
      HPV_pos == 1 ~ "Yes"
    ), 
    
    # sociodemographic characteristics
    AGE = as.numeric(AGE),
    M.STATUS = case_when(
      M.STATUS == "MARRIED" ~ "Married",
      M.STATUS == "SEPARATED" ~ "Married", 
      M.STATUS == "SINGLE" ~ "Single",
      TRUE ~ NA_character_
    ), 
    EDU.LEVEL = case_when(
      EDU.LEVEL == "PRIMARY" ~ "Primary", 
      EDU.LEVEL == "SECONDARY" ~ "Secondary", 
      EDU.LEVEL == "TERTIARY" ~ "Tertiary",
      TRUE ~ NA_character_
    ), 
    hpv_coinfect = ifelse(Number_of_HPV_types == 0, NA, Number_of_HPV_types), 
    
    # HIV and STI infections
    HIV = case_when(
      HIV == 0 ~ "No", 
      HIV == 1 ~ "Yes", 
    ), 
    HIV_code1 = case_when(
      HIV_code1 == "HIV-" ~ NA_character_, 
      HIV_code1 == "HIV+" ~ "Yes", 
      HIV_code1 == "ART-" ~ "No"
    ),
    HIV_code1 = factor(HIV_code1), 
    HIV_code1 = fct_relevel(HIV_code1, "Yes"), 
    HIV_code1_all = case_when(
      HIV_code1 == "ART-" ~ "No ART",
      HIV_code1 == "HIV-" & PREP == "Y" ~ "PREP",
      HIV_code1 == "HIV-" & PEP == "Y" ~ "PEP",
      HIV_code1 == "HIV+" ~ "ART"
    ), 
    CD4C_cat = case_when(
      CD4C <= 500 ~ "<=500", 
      CD4C > 500 ~ "500+"
    ), 
    VL = case_when(
      VL == 40 ~ "Yes",
      VL != 40 ~ "No", 
      TRUE ~ NA_character_
    ), 
    STI.RESULT = case_when(
      STI.RESULT == "STI" ~ "Yes", 
      STI.RESULT == "NO_STI" ~ "No", 
      TRUE ~ NA_character_
    ), 
    STI_RESULT2 = case_when(
      STI_RESULT2 == "NA" ~ NA_character_,
      STI_RESULT2 == "NO_STI" ~ "No", 
      STI_RESULT2 == "OTHER" ~ "Other", 
      STI_RESULT2 == "G/C" ~ "Gonorrhea", 
      TRUE ~ NA_character_
    ), 
    STI_RESULT2 = factor(STI_RESULT2), 
    STI_RESULT2 = fct_relevel(STI_RESULT2, "Gonorrhea")
  ) 


collapse_data <- function(data, vars, name) {
  data %>%
    select(PID, {{ vars }}) %>%
    mutate(across({{ vars }}, ~ as.numeric(.x))) %>%
    pivot_longer(-PID) %>%
    filter(!is.na(value)) %>%
    group_by(PID) %>% 
    summarise(n = sum(value, na.rm = TRUE)) %>%
    {
      names(.)[2] <- name
      .
    }
}

# sexual practices and hygiene 
reg <- collapse_data(
  data.sd,
  c(REGULAR_partner, CASUAL_partner, P.P), 
  "all_partner"
)  
oral <- collapse_data(
  data.sd, 
  c(ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay), 
  "oral_partner"
)
inset <- collapse_data(
  data.sd, 
  c(INSERT_sex_reg, INSERT_sex_cas, INSERT_sex_pay), 
  "insert_partner"
)
recep <- collapse_data(
  data.sd, 
  c(RECEPT_sex_reg, RECEPT_sex_cas, RECEPT_sex_pay), 
  "recep_partner"
)
cond <- collapse_data(
  data.sd,
  c(CONDOM_reg, CONDOM_cas, CONDOM_pay), 
  "condom_partner"
)

hpv <- hpv %>% 
  left_join(reg, by = "PID") %>% 
  left_join(oral, by = "PID") %>% 
  left_join(inset, by = "PID") %>% 
  left_join(recep, by = "PID") %>% 
  left_join(cond, by = "PID") %>% 
  mutate(
    MSM.STATUS.x = case_when(
      MSM.STATUS.x == "MSM" ~ "Yes", 
      MSM.STATUS.x == "NON-MSM" ~ "No", 
      TRUE ~ NA_character_
    ), 
    S.WORKER = case_when(
      S.WORKER == "Y" ~ "Yes", 
      S.WORKER == "N" ~ "No", 
      TRUE ~ NA_character_ 
    ), 
    D.PROST = as.numeric(D.PROST), 
    REGULAR_partner = as.numeric(REGULAR_partner),
    CASUAL_partner = as.numeric(CASUAL_partner),
    P.P = as.numeric(P.P),
    ORAL_sex_reg = as.numeric(ORAL_sex_reg),
    ORAL_sex_cas = as.numeric(ORAL_sex_cas),
    ORAL_sex_pay = as.numeric(ORAL_sex_pay),
    INSERT_sex_reg = as.numeric(INSERT_sex_reg),
    INSERT_sex_cas = as.numeric(INSERT_sex_cas),
    INSERT_sex_pay = as.numeric(INSERT_sex_pay),
    RECEPT_sex_reg = as.numeric(RECEPT_sex_reg),
    RECEPT_sex_cas = as.numeric(RECEPT_sex_cas),
    RECEPT_sex_pay = as.numeric(RECEPT_sex_pay),
    CONDOM_reg = as.numeric(CONDOM_reg),
    CONDOM_cas = as.numeric(CONDOM_cas),
    CONDOM_pay = as.numeric(CONDOM_pay),
    DOUCHING = case_when(
      DOUCHING == "Y" ~ "Yes", 
      DOUCHING == "N" ~ "No", 
      TRUE ~ NA_character_
    ), 
    AGE = set_label(AGE, "Age"), 
    M.STATUS = set_label(M.STATUS, "Marital status"), 
    EDU.LEVEL = set_label(EDU.LEVEL, "Education level"), 
    hpv_coinfect = set_label(hpv_coinfect, "Multiple HPV infection"), 
    HIV = set_label(HIV, "HIV status"), 
    HIV_code1 = set_label(HIV_code1, "ART status"), 
    CD4C_cat = set_label(CD4C_cat, "CD4 count category"), 
    CD4pc = set_label(CD4pc, "CD4 percentage"), 
    VL = set_label(VL, "Undetectable viral load"), 
    STI_RESULT2 = set_label(STI_RESULT2, "STI infection"),
    STI.RESULT = set_label(STI.RESULT, "STI infection"),
    
    MSM.STATUS.x = set_label(MSM.STATUS.x, "MSM"), 
    S.WORKER = set_label(S.WORKER, "Sex worker"), 
    D.PROST = set_label(D.PROST, "Duration of sex work"), 
    all_partner = set_label(all_partner, "Partners"), 
    REGULAR_partner = set_label(REGULAR_partner, "Regular"),
    CASUAL_partner = set_label(CASUAL_partner, "Casual"), 
    P.P = set_label(P.P, "Transactional"),
    oral_partner = set_label(oral_partner, "Oral sex with partners"), 
    ORAL_sex_reg = set_label(ORAL_sex_reg, "Regular"), 
    ORAL_sex_cas = set_label(ORAL_sex_cas, "Casual"), 
    ORAL_sex_pay = set_label(ORAL_sex_pay, "Transactional"), 
    insert_partner = set_label(insert_partner, "Insertive sex with partners"), 
    INSERT_sex_reg = set_label(INSERT_sex_reg, "Regular"), 
    INSERT_sex_cas = set_label(INSERT_sex_cas, "Casual"), 
    INSERT_sex_pay = set_label(INSERT_sex_pay, "Transactional"), 
    recep_partner = set_label(recep_partner, "Receptive sex with partners"), 
    RECEPT_sex_reg = set_label(RECEPT_sex_reg, "Regular"), 
    RECEPT_sex_cas = set_label(RECEPT_sex_reg, "Casual"), 
    RECEPT_sex_pay = set_label(RECEPT_sex_reg, "Transactional"), 
    condom_partner = set_label(condom_partner, "Sex with condom use"), 
    CONDOM_reg = set_label(CONDOM_reg, "Regular"), 
    CONDOM_cas = set_label(CONDOM_cas, "Casual"), 
    CONDOM_pay = set_label(CONDOM_pay, "Transactional"), 
    DOUCHING = set_label(DOUCHING, "Douching practice")
) 



# Table 1: Socio-demographic characteristics ------------------------------

t1 <- hpv %>% 
  tbl_summary(
    by = HPV_pos_cat, 
    include = c(AGE, M.STATUS, EDU.LEVEL, hpv_coinfect), 
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p})"), 
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ 1
    ),
    missing = "no"
  ) %>% 
  add_overall() %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HPV Infection**") %>%
  as_flex_table() %>% 
  compose(i = (10:14) + 0, j = 3, as_paragraph(as_chunk("-"))) %>% 
  compose(i = 9 + 0, j = 5, as_paragraph(as_chunk("-"))) %>% 
  set_caption("Clients' socio-demographic characteristics by HPV infection")
t1

# Table 2: HIV infection and related information
t2 <- hpv %>% 
  tbl_summary(
    by = HPV_pos_cat, 
    include = c(HIV, HIV_code1, CD4C_cat, CD4pc, VL, STI_RESULT2), 
    statistic = list(c(CD4pc) ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ 1
    ),
    missing = "no"
  ) %>% 
  add_overall() %>%
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HPV Infection**") %>%
  as_flex_table() %>% 
  set_caption("Clients' HIV status and related information by HPV status")
t2


## create table 3
t3 <- hpv %>%
  select(
    HPV_pos_cat, MSM.STATUS.x, S.WORKER, D.PROST,
    all_partner, REGULAR_partner, CASUAL_partner, P.P, 
    oral_partner, ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay, 
    insert_partner, INSERT_sex_reg, INSERT_sex_cas,
    INSERT_sex_pay, recep_partner, RECEPT_sex_reg, 
    RECEPT_sex_cas, RECEPT_sex_pay, condom_partner,
    CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING
  ) %>% 
  tbl_summary(
    by = HPV_pos_cat,
    type = list(c(
      all_partner, REGULAR_partner, CASUAL_partner, P.P, 
      oral_partner, ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay, 
      insert_partner, INSERT_sex_reg, INSERT_sex_cas,
      INSERT_sex_pay, recep_partner, RECEPT_sex_reg, 
      RECEPT_sex_cas, RECEPT_sex_pay, condom_partner,
      CONDOM_reg, CONDOM_cas, CONDOM_pay
    ) ~ "continuous"),
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ 1
    ),
    missing = "no"
  ) %>% 
  add_overall() %>%
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HPV Infection**") %>%
  as_flex_table() %>% 
  bg(i = c(1, 5, 9, 13, 17, 21) + 3, bg = "lightgrey") %>%
  set_caption("Clients' sexual behaviors and hygiene practices by HPV status")
t3 




# HIV genotype distribution -----------------------------------------------

genotype <- nml %>% 
  left_join(
    hpv %>% 
      select(PID, AGE, MSM.STATUS.x, S.WORKER, HIV, 
             B.YEAR,	AGE, AGE.F.SEX,	M.STATUS,	EDU.LEVEL, 
             CD4pc, CD4C_cat, VL, HIV_code1, PREP, PEP,
             STI_TEST, STI_RESULT2,	STI.RESULT, 
             REGULAR_partner, CASUAL_partner, P.P, ORAL_sex_reg, 
             ORAL_sex_cas, ORAL_sex_pay, INSERT_sex_reg, INSERT_sex_cas,
             INSERT_sex_pay, RECEPT_sex_reg, RECEPT_sex_cas, RECEPT_sex_pay,
             CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING, DRUG.HIST, 
             hpv_coinfect), 
    by = "PID"
  ) %>% 
  mutate(
    age_cat = ifelse(AGE <= 25, "<=25", "25+"), 
    hr_hpv_all = ifelse(hr_hpv_all == "", NA, hr_hpv_all), 
    lr_hpv_all = ifelse(lr_hpv_all == "", NA, lr_hpv_all)
  )

ts1 <- genotype %>% 
  tbl_summary(
    by = age_cat, 
    include = num_hpv_all, 
    statistic = list(all_categorical() ~ "{n} ({p})"), 
    label = list(num_hpv_all ~ "HPV Type")
  ) %>% 
  add_overall() %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Age**") %>%
  as_flex_table()  %>%
  set_caption("HPV genotypes by age categories")
ts1

t4 <- genotype %>% 
  tbl_summary(
    by = age_cat, 
    include = c(hpv_all, hpv_16, hpv_18, hr_hpv_all, lr_hpv_all), 
    statistic = list(all_categorical() ~ "{n} ({p})"), 
    label = list(
      hpv_all ~ "HPV Type", 
      hpv_16 ~ "HPV 16", 
      hpv_18 ~ "HPV 18", 
      hr_hpv_all ~ "High risk HPV Type", 
      lr_hpv_all ~ "Low risk HPV Type"
    )
  ) %>% 
  add_overall() %>% 
  # add_p() %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Age**") %>%
  as_flex_table()  %>%
  set_caption("Different HPV genotypes by age categories")
t4

t5 <- genotype %>% 
  tbl_summary(
    by = HIV, 
    include = c(hpv_all, hpv_16, hpv_18, hr_hpv_all, lr_hpv_all), 
    statistic = list(all_categorical() ~ "{n} ({p})"), 
    label = list(
      hpv_all ~ "HPV Type", 
      hpv_16 ~ "HPV 16", 
      hpv_18 ~ "HPV 18", 
      hr_hpv_all ~ "High risk HPV Type", 
      lr_hpv_all ~ "Low risk HPV Type"
    ), 
    missing = "no"
  ) %>% 
  add_overall() %>% 
  # add_p() %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HIV Status**") %>%
  as_flex_table() 
t5 

create_barplot <- function(data, var, title) {
  data %>% 
    group_by(hpv_all, {{ var }}) %>% 
    summarize(n = n()) %>%
    filter(!is.na({{ var }})) %>% 
    ggplot(aes({{ var }}, n, fill = hpv_all)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "lightgrey")) + 
    ylab("Number") + 
    xlab(title) +
    theme_classic()
}

p_hpv_coinfect <- create_barplot(genotype, hpv_coinfect, "Multiple HPV infection") 
p_age <- create_barplot(genotype, age_cat, "Age category")
p_marital <- create_barplot(genotype, M.STATUS, "Marital Status")
p_edu <- create_barplot(genotype, EDU.LEVEL, "Education Level")
p_hiv <- create_barplot(genotype, HIV, "HIV Status")
p_cd4c <- create_barplot(genotype, CD4C_cat, "CD4 Count")
p_vl <- create_barplot(genotype, VL, "Undetectable Viral Load")
p_sti <- create_barplot(genotype, STI_RESULT2, "STI infection")
p_msm <- create_barplot(genotype, MSM.STATUS.x, "MSM")
p_sexwork <- create_barplot(genotype, S.WORKER, "Sex Worker")

library(patchwork)
plot <- p_hpv_coinfect + 
  p_age + 
  p_marital + 
  p_edu + 
  p_hiv + 
  p_cd4c + 
  p_vl + 
  p_sti + 
  p_msm + 
  p_sexwork + 
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  guides(fill=guide_legend(title="HPV Type")) & 
  theme(legend.position = "bottom", plot.tag = element_text(size = 8))  
plot




# Regression analyses -----------------------------------------------------


hpv_reg <- hpv %>% 
  mutate(age_cat = ifelse(AGE <= 25, "<=25", "25+")) %>% 
  select(HPV_pos, age_cat, M.STATUS, EDU.LEVEL, 
           HIV, CD4pc, STI.RESULT, all_partner, 
           insert_partner, recep_partner, condom_partner) %>% 
  na.omit()




run_glm <- function(x_name) {
  f <- sprintf("HPV_pos ~ %s", paste(x_name, collapse = " + "))
  m <- glm(formula(f), data = hpv_reg, family = binomial)
  tbl_regression(m, exponentiate = TRUE, 
                 pvalue_fun = purrr::partial(style_pvalue, digits = 3))
} 

t_age <- run_glm("age_cat")
t_marital <- run_glm("M.STATUS")
t_edu <- run_glm("EDU.LEVEL")
t_HIV <- run_glm("HIV")
t_cd4pc <- run_glm("CD4pc")
t_sti <- run_glm("STI.RESULT")
t_allp <- run_glm("all_partner")
t_insert <- run_glm("insert_partner")
t_recep <- run_glm("recep_partner")
t_condom <- run_glm("condom_partner")

# combine unadjusted analyses
t_unadj <- tbl_stack(list(t_age, t_marital, t_edu, t_HIV, t_cd4pc, 
                    t_sti, t_allp, t_insert, t_recep, t_condom))  
t_unadj

# run multivariable logistic regression 
m <- glm(HPV_pos ~ age_cat + M.STATUS + EDU.LEVEL + 
           HIV + CD4pc + STI.RESULT + all_partner + 
           insert_partner + recep_partner + condom_partner,
         hpv_reg, family = binomial)
# run backward step-wise regression 
step(m, direction = "backward")

# final model 
mf <- glm(formula = HPV_pos ~ M.STATUS + HIV + insert_partner + recep_partner, 
    family = binomial, data = hpv_reg)

t_adj <- tbl_regression(mf, exponentiate = TRUE, 
                        pvalue_fun = purrr::partial(style_pvalue, digits = 3))
t_adj

t_regression <- tbl_merge(list(t_unadj, t_adj), 
          tab_spanner = c("**Unadjusted**", "**Adjusted**")) %>% 
  as_flex_table() %>%
  set_caption("Association between risk factors and HPV infection, N = 70")
  
t_regression











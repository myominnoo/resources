

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#                 HPV Microbiome Data Analysis
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Author: Myo Minn Oo
# Start Date: 30-06-2022
# Revised Date: 09-08-2022
# Revised Date: 24-08-2022
# Revised Date: 06-09-2022

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
# 81  11  16  18  26  33  35  42  45  51  52  54  56  58  59   6  66  67  68  
# 1   8  10   2   1   2   7   4   7   3   1   1   4   6   6   9   2   5   1   
# 69  70  72  73  81  83  84  90  97 
# 3   5   3   2   3   3   3   1   1

hpv_all <- c(81, 11, 16, 18, 26, 33, 35, 42, 45, 51, 52, 54, 56, 58, 59, 6, 66, 67, 68, 
             69, 70, 72, 73, 81, 83, 84, 90, 97)
# set high risk and low risk hpv types
# hr_hpv_types <- c(31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 
#                   26, 53, 73, 82, 66, 68)
# lr_hpv_types <- c(6, 11, 30, 67, 70)
## ref: Bouassa et al 2018 
# HR: 16, -18, -31, -33, -35, -39, -45, -51, -52, -56, -58, -59, and -68
# LR: -6, -11, -40, -42, -43, -44,-53, -54 and -70
# possibly carcinogenic: -26, -61, -66, -69, -73 and -82
hr_hpv_types <- c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68)
cancer_hpv <- c(26, 61, 69, 73, 82)
# lr_hpv_types <- c(6, 11, 40, 42, 43, 44, 53, 54, 70)
lr_hpv_types <- sort(setdiff(hpv_all, c(hr_hpv_types, cancer_hpv)))

# The 4-valent Gardasil-41 vaccine (Merck & Co. Inc., New Jersey, USA) is 
# effective against HPV genotypes 6, 11, 16 and 18;
val4_hpv <- c(6, 11, 16, 18)

# The 9-valent Gardasil-91 vaccine iseffective against HPV genotypes 
# 6, 11, 16, 18, 31, 33, 45, 52 and 58.
val9_hpv <- c(6, 11, 16, 18, 31, 33, 45, 52, 58)
hpv_16_18 <- c(16, 18)

genotype_all <- data.sd %>%
  select(PID, NMLTypingResult) %>% 
  filter(!grepl("nega|incon", NMLTypingResult)) %>%
  {
    nmax <- max(stringr::str_count(.$NMLTypingResult, ","), na.rm = TRUE) + 1
    separate(., NMLTypingResult, paste0("hpv.", seq_len(nmax)), 
             sep = ",", fill = "right")
  } %>% 
  pivot_longer(-PID) %>%
  select(-name) %>% 
  filter(!is.na(value)) 

nml <- genotype_all %>% 
  mutate(value = as.numeric(value), 
         name = paste0("hpv_", value)) %>% 
  arrange(value) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(num_hpv = rowSums(!across(contains("hpv"), is.na)),
         hr_hpv = rowSums(!across(num_range("hpv_", hr_hpv_types), is.na)),
         lr_hpv = rowSums(!across(num_range("hpv_", lr_hpv_types), is.na)),
         cancer_hpv = rowSums(!across(num_range("hpv_", cancer_hpv), is.na)),
         val4_hpv = rowSums(!across(num_range("hpv_", val4_hpv), is.na)),
         val9_hpv = rowSums(!across(num_range("hpv_", val9_hpv), is.na)),
         hpv1618 = rowSums(!across(num_range("hpv_", hpv_16_18), is.na))) %>%

  unite(num_hpv_all, contains("hpv_"), sep = "+", na.rm = TRUE, remove = FALSE) %>% 
  unite(hr_hpv_all, num_range("hpv_", hr_hpv_types), sep = "+",
        na.rm = TRUE, remove = FALSE) %>% 
  unite(lr_hpv_all, num_range("hpv_", lr_hpv_types), sep = "+",
        na.rm = TRUE, remove = FALSE) %>% 
  mutate(hpv_16 = ifelse(is.na(hpv_16), 0, 1), 
         hpv_18 = ifelse(is.na(hpv_18), 0, 1), 
         hpv_all = case_when(
           hpv_16 > 0 | hpv_18 > 0 ~ "HPV 16+18", 
           hr_hpv > 0 ~ "HR HPV", 
           lr_hpv > 0 ~ "LR HPV", 
           TRUE ~ "Other HPV"
         )) %>% 
  mutate(
    any_hpv = ifelse(num_hpv > 0, 1, 0),
    multi_hpv = ifelse(num_hpv > 1, 1, 0), 
    any_hr_hpv = ifelse(hr_hpv > 0, 1, 0), 
    any_lr_hpv = ifelse(lr_hpv > 0, 1, 0), 
    any_cancer_hpv = ifelse(cancer_hpv > 0, 1, 0), 
    multi_hr_hpv = ifelse(hr_hpv > 1, 1, 0),
    hpv1618 = ifelse(hpv1618 > 0, 1, 0), 
    any_val4 = ifelse(val4_hpv > 0, 1, 0), 
    multi_val4 = ifelse(val4_hpv > 1, 1, 0), 
    any_val9 = ifelse(val9_hpv > 0, 1, 0), 
    multi_val9 = ifelse(val9_hpv > 1, 1, 0)
  ) %>% 
  select(PID, num_hpv_all, lr_hpv_all, hr_hpv_all, hpv_16, hpv_18, hpv1618, 
         num_hpv:multi_val9) 



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
    age_cat = case_when(
      AGE <= 20 ~ 1, 
      AGE > 20 & AGE <=30 ~ 2, 
      AGE > 30 ~ 3,
      TRUE ~ NA_real_
    ), 
    age_cat = factor(age_cat, labels = c("<=20", "21-30", "31+")),
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
      HIV == 0 ~ "Negative", 
      HIV == 1 ~ "Positive", 
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
      STI_RESULT2 == "NA" ~ NA_real_,
      STI_RESULT2 == "NO_STI" ~ 1, 
      STI_RESULT2 == "OTHER" ~ 3, 
      STI_RESULT2 == "G/C" ~ 2, 
      TRUE ~ NA_real_
    ), 
    STI_RESULT2 = factor(STI_RESULT2, labels = c("No", "Gonorrhea", "Other")), 
    num_sex_partner = as.numeric(REGULAR_partner) + as.numeric(CASUAL_partner) + 
      as.numeric(P.P)
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
    # MSM.STATUS.x = case_when(
    #   MSM.STATUS.x == "MSM" ~ "Yes",
    #   MSM.STATUS.x == "NON-MSM" ~ "No",
    #   TRUE ~ NA_character_
    # ),
    MSM.STATUS.x = factor(MSM.STATUS.x, levels = c("NON-MSM", "MSM")), 
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
    across(c(
      all_partner, REGULAR_partner, CASUAL_partner, P.P, 
      oral_partner, ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay, 
      insert_partner, INSERT_sex_reg, INSERT_sex_cas,
      INSERT_sex_pay, recep_partner, RECEPT_sex_reg, 
      RECEPT_sex_cas, RECEPT_sex_pay, condom_partner,
      CONDOM_reg, CONDOM_cas, CONDOM_pay
    ), ~ case_when(
      is.na(.x) | .x == 0 ~ "No", 
      TRUE ~ "Yes"
    )), 
    
    AGE = set_label(AGE, "Age"), 
    age_cat = set_label(age_cat, "Age category"), 
    M.STATUS = set_label(M.STATUS, "Marital status"), 
    EDU.LEVEL = set_label(EDU.LEVEL, "Education level"), 
    hpv_coinfect = set_label(hpv_coinfect, "Multiple HPV infection"), 
    HIV = set_label(HIV, "HIV status"), 
    HIV_code1 = set_label(HIV_code1, "ART status"), 
    CD4C_cat = set_label(CD4C_cat, "CD4 count category"), 
    CD4pc = set_label(CD4pc, "CD4 %"), 
    VL = set_label(VL, "Undetectable viral load"), 
    STI_RESULT2 = set_label(STI_RESULT2, "STI infection"),
    STI.RESULT = set_label(STI.RESULT, "STI infection"),
    
    MSM.STATUS.x = set_label(MSM.STATUS.x, "MSM"), 
    S.WORKER = set_label(S.WORKER, "Sex worker"), 
    D.PROST = set_label(D.PROST, "Duration of sex work"), 
    num_sex_partner = set_label(num_sex_partner, "Number of sex partners"),
    all_partner = set_label(all_partner, "All sex partners"), 
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
    DOUCHING = set_label(DOUCHING, "Douching practice"), 
    
    nocondom_recept_reg = ifelse(CONDOM_reg == "No" & RECEPT_sex_reg == "Yes", 1, 0), 
    nocondom_recept_cas = ifelse(CONDOM_cas == "No" & RECEPT_sex_cas == "Yes", 1, 0), 
    nocondom_recept_pay = ifelse(CONDOM_pay == "No" & RECEPT_sex_pay == "Yes", 1, 0), 
    nocondom_insert_reg = ifelse(CONDOM_reg == "No" & INSERT_sex_reg == "Yes", 1, 0), 
    nocondom_insert_cas = ifelse(CONDOM_cas == "No" & INSERT_sex_cas == "Yes", 1, 0), 
    nocondom_insert_pay = ifelse(CONDOM_pay == "No" & INSERT_sex_pay == "Yes", 1, 0), 
    nocondom_recept = ifelse(condom_partner == "No" & recep_partner == "Yes", 1, 0), 
    nocondom_insert = ifelse(condom_partner == "No" & insert_partner == "Yes", 1, 0), 
    
    nocondom_recept = set_label(nocondom_recept, "Condomless receptive sex"), 
    nocondom_recept_reg = set_label(nocondom_recept_reg, "Regular"), 
    nocondom_recept_cas = set_label(nocondom_recept_cas, "Casual"), 
    nocondom_recept_pay = set_label(nocondom_recept_pay, "Transactional"), 
    
    nocondom_insert = set_label(nocondom_insert, "Condomless insertive sex"), 
    nocondom_insert_reg = set_label(nocondom_insert_reg, "Regular"),
    nocondom_insert_cas = set_label(nocondom_insert_cas, "Casual"), 
    nocondom_insert_pay = set_label(nocondom_insert_pay, "Transactional")
) 


summary(hpv$AGE)

hpv %>% 
  group_by(MSM.STATUS.x) %>% 
  summarize(mean = mean(AGE, na.rm = TRUE), 
            sd = sd(AGE, TRUE))

hpv %>% 
  
  tbl_summary(
    by = MSM.STATUS.x,
    include = age_cat, 
    digits = all_categorical() ~ c(0, 1)
  )

# Table 1: Socio-demographic characteristics ------------------------------

t1_total <- hpv %>% 
  tbl_summary(
    include = c(AGE, age_cat, M.STATUS, EDU.LEVEL), 
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p})"), 
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "column"
  ) %>% 
  modify_header(update = stat_0 ~ "**n(%), N = {N}**")

t1_msm <- hpv %>% 
  tbl_summary(
    by = MSM.STATUS.x, 
    include = c(AGE, age_cat, M.STATUS, EDU.LEVEL), 
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p})"), 
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**MSM Status**") 



t1_hiv <- hpv %>% 
  tbl_summary(
    by = HIV, 
    include = c(AGE, age_cat, M.STATUS, EDU.LEVEL), 
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p})"), 
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) %>% 
  # adding spanning header
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HIV Status**") 

t1 <- tbl_merge(tbls = list(t1_total, t1_msm, t1_hiv), 
          tab_spanner = c("Total", "MSM", "HIV Status")) %>%
  as_flex_table() %>% 
  compose(i = c(1, 3:5), j = 4, as_paragraph(as_chunk("-"))) %>% 
  compose(i = 2, j = 5, as_paragraph(as_chunk("-"))) %>% 
  set_caption("socio-demographic characteristics by sexual orientation and HIV status")
t1


# Table 2: HIV infection and related information

t2_total <- hpv %>% 
  tbl_summary(
    include = c(HIV, HIV_code1, CD4C_cat, CD4pc, VL, STI_RESULT2), 
    statistic = list(c(CD4pc) ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "column"
  ) %>% 
  modify_header(update = stat_0 ~ "**n(%), N = {N}**")


t2_msm <- hpv %>% 
  tbl_summary(
    by = MSM.STATUS.x, 
    include = c(HIV, HIV_code1, CD4C_cat, CD4pc, VL, STI_RESULT2), 
    statistic = list(c(CD4pc) ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) 

t2_hiv <- hpv %>% 
  tbl_summary(
    by = HIV, 
    include = c(HIV, HIV_code1, CD4C_cat, CD4pc, VL, STI_RESULT2), 
    statistic = list(c(CD4pc) ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) 

t2 <- tbl_merge(tbls = list(t2_total, t2_msm, t2_hiv), 
          tab_spanner = c("Total", "MSM", "HIV Status")) %>%
  as_flex_table() %>% 
  compose(i = 1, j = 8, as_paragraph(as_chunk("-"))) %>% 
  compose(i = 2:3, j = 6:7, as_paragraph(as_chunk("-"))) %>% 
  compose(i = 5, j = 8, as_paragraph(as_chunk("-"))) %>% 
  compose(i = c(4, 6:9), j = c(6, 8), as_paragraph(as_chunk("-"))) %>% 
  set_caption("HIV status and related information by sexual orientation and HIV status")
t2


# Table 3 -----------------------------------------------------------------

t3_total <- hpv %>%
  tbl_summary(
    include = c(
      MSM.STATUS.x, S.WORKER, D.PROST, num_sex_partner,
      all_partner, REGULAR_partner, CASUAL_partner, P.P, 
      oral_partner, ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay, 
      insert_partner, INSERT_sex_reg, INSERT_sex_cas,
      INSERT_sex_pay, recep_partner, RECEPT_sex_reg, 
      RECEPT_sex_cas, RECEPT_sex_pay, condom_partner,
      CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING
    ), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "column"
  ) %>% 
  modify_header(update = stat_0 ~ "**n(%), N = {N}**")

t3_msm <- hpv %>%
  tbl_summary(
    by = MSM.STATUS.x, 
    include = c(
      MSM.STATUS.x, S.WORKER, D.PROST, num_sex_partner,
      all_partner, REGULAR_partner, CASUAL_partner, P.P, 
      oral_partner, ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay, 
      insert_partner, INSERT_sex_reg, INSERT_sex_cas,
      INSERT_sex_pay, recep_partner, RECEPT_sex_reg, 
      RECEPT_sex_cas, RECEPT_sex_pay, condom_partner,
      CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING
    ), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) 

t3_hiv <- hpv %>%
  tbl_summary(
    by = HIV, 
    include = c(
      MSM.STATUS.x, S.WORKER, D.PROST, num_sex_partner,
      all_partner, REGULAR_partner, CASUAL_partner, P.P, 
      oral_partner, ORAL_sex_reg, ORAL_sex_cas, ORAL_sex_pay, 
      insert_partner, INSERT_sex_reg, INSERT_sex_cas,
      INSERT_sex_pay, recep_partner, RECEPT_sex_reg, 
      RECEPT_sex_cas, RECEPT_sex_pay, condom_partner,
      CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING
    ), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) 


t3 <- tbl_merge(tbls = list(t3_total, t3_msm, t3_hiv), 
          tab_spanner = c("Total", "MSM", "HIV Status")) %>%
  as_flex_table() %>% 
  bg(i = c(1, 5, 9, 13, 17, 21) + 6, bg = "lightgrey") %>%
  set_caption("HIV sexual behaviors and hygiene practices by sexual orientation and HIV status")
t3 



# HIV genotype distribution -----------------------------------------------

# Figure 1: distribution of anal hpv types by 4- and 9- valent vaccines
genotype_all <- genotype_all %>% 
  left_join(
    hpv %>% 
      select(PID, HIV, MSM.STATUS.x), 
    by = "PID"
  ) %>% 
  mutate(value = as.numeric(value), 
         type = factor(value), 
         `Valent Vaccine` = ifelse(value %in% c(val4_hpv, val9_hpv), 
                       "4- or 9-valent vaccine type", "Non-HPV vaccine type"), 
         hpv = case_when(
           value %in% hr_hpv_types ~ "HR-HPV", 
           value %in% lr_hpv_types ~ "LR-HPV", 
           value %in% cancer_hpv ~ "Possibly oncogenic", 
           TRUE ~ NA_character_
         )) 
  
fig1a <- genotype_all %>% 
  group_by(hpv, `Valent Vaccine`, type) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(type, n, fill = `Valent Vaccine`)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  scale_y_continuous(name = "Number", breaks = seq(0, 12, 2)) + 
  facet_wrap(~ hpv, scales = "free_x") + 
  theme_classic()

fig1b <- genotype_all %>% 
  group_by(hpv, HIV, type) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(type, n, fill = HIV)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  scale_y_continuous(name = "Number", breaks = seq(0, 12, 2)) + 
  facet_wrap(~ hpv, scales = "free_x") + 
  theme_classic()

fig1c <- genotype_all %>% 
  group_by(hpv, MSM.STATUS.x, type) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(type, n, fill = MSM.STATUS.x)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  scale_y_continuous(name = "Number", breaks = seq(0, 12, 2)) + 
  facet_wrap(~ hpv, scales = "free_x") + 
  theme_classic()

library(patchwork)

fig1 <- fig1a / fig1c / fig1b + 
  plot_annotation(tag_levels = "A")



# Table 4: ----------------------------------------------------------------

genotype <- hpv %>% 
  select(PID, AGE, MSM.STATUS.x, S.WORKER, HIV, 
         B.YEAR,	AGE, AGE.F.SEX,	M.STATUS,	EDU.LEVEL, 
         CD4pc, CD4C_cat, VL, HIV_code1, PREP, PEP,
         STI_TEST, STI_RESULT2,	STI.RESULT, 
         num_sex_partner, all_partner, 
         insert_partner, recep_partner, condom_partner, 
         REGULAR_partner, CASUAL_partner, P.P, ORAL_sex_reg, 
         ORAL_sex_cas, ORAL_sex_pay, INSERT_sex_reg, INSERT_sex_cas,
         INSERT_sex_pay, RECEPT_sex_reg, RECEPT_sex_cas, RECEPT_sex_pay,
         CONDOM_reg, CONDOM_cas, CONDOM_pay, DOUCHING, DRUG.HIST, 
         hpv_coinfect) %>%
  left_join(nml, by = "PID") %>%
  mutate(
    age_cat = ifelse(AGE <= 25, "<=25", "25+"), 
    age_cat = set_label(age_cat, "Age category"), 
    hr_hpv_all = ifelse(hr_hpv_all == "", NA, hr_hpv_all), 
    lr_hpv_all = ifelse(lr_hpv_all == "", NA, lr_hpv_all), 
    any_hpv = ifelse(is.na(any_hpv), 0, any_hpv), 
    across(hpv_coinfect:multi_hr_hpv, ~ ifelse(is.na(.x), 0, .x)),
    
    any_hpv = set_label(any_hpv, "Any HPV"), 
    multi_hpv = set_label(multi_hpv, "Multiple type of any HPV"), 
    any_lr_hpv = set_label(any_lr_hpv, "LR-HPV"), 
    any_hr_hpv = set_label(any_hr_hpv, "HR-HPV"), 
    multi_hr_hpv = set_label(multi_hr_hpv, "Multiple type of HR-HPV"), 
    hpv_16 = set_label(hpv_16, "HPV-16"), 
    hpv_18 = set_label(hpv_18, "HPV-18"), 
    hpv1618 = set_label(hpv1618, "HPV-16 and HPV-18"), 
    any_val4 = set_label(any_val4, "Any 4-valent vaccine type"), 
    multi_val4 = set_label(multi_val4, "Multiple 4-valent vaccine type"), 
    any_val9 = set_label(any_val9, "Any 9-valent vaccine type"), 
    multi_val9 = set_label(multi_val9, "Multiple 9-valent vaccine type"), 
    
    num_sex_partner = case_when(
      num_sex_partner > 0 & num_sex_partner <= 5 ~ "1-5", 
      num_sex_partner > 5 ~ "5+", 
      TRUE ~ "0"
    ), 
    num_sex_partner = set_label(num_sex_partner, "Number of sex partners")
  ) 

t4_total <- genotype %>% 
  filter(any_hpv == 1) %>% 
  tbl_summary(
    include = c(multi_hpv, any_lr_hpv, 
                any_hr_hpv, multi_hr_hpv, hpv_16, hpv_18, hpv1618, 
                any_val4, multi_val4, any_val9, multi_val9), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "column"
  ) %>% 
  modify_header(update = stat_0 ~ "**n(%), N = {N}**")

t4_msm <- genotype %>%
  filter(any_hpv == 1) %>% 
  tbl_summary(
    by = MSM.STATUS.x, 
    include = c(multi_hpv, any_lr_hpv, 
                any_hr_hpv, multi_hr_hpv, hpv_16, hpv_18, hpv1618, 
                any_val4, multi_val4, any_val9, multi_val9), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) 

t4_hiv <- genotype %>%
  filter(any_hpv == 1) %>% 
  tbl_summary(
    by = HIV, 
    include = c(multi_hpv, any_lr_hpv, 
                any_hr_hpv, multi_hr_hpv, hpv_16, hpv_18, hpv1618, 
                any_val4, multi_val4, any_val9, multi_val9), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "row"
  ) %>% 
  # rounding p-values to 3 decimal places
  add_p(pvalue_fun = function(x) style_number(x, digits = 3)) 

t4 <- tbl_merge(tbls = list(t4_total, t4_msm, t4_hiv), 
                tab_spanner = c("Total", "MSM", "HIV Status")) %>%
  as_flex_table()  %>%
  set_caption("Anal HPV genotype distributions by sexual orientation and HIV status")
t4



# Regression analyses -----------------------------------------------------

glm_uni <- function(y, x, data, stats = "**cOR (95% CI)**") {
  f <- formula(paste(y, " ~ ", x))
  glm(f, data = data, family = binomial) %>%
    tbl_regression(
      exponentiate = TRUE, 
      estimate_fun = purrr::partial(style_ratio, digits = 1),
      pvalue_fun = purrr::partial(style_pvalue, digits = 3)
    ) %>% 
    # merge OR and CI into single column
    modify_table_styling(
      column = estimate,
      rows = !is.na(estimate),
      cols_merge_pattern = "{estimate} ({conf.low}, {conf.high})"
    ) %>%
    modify_header(estimate ~ stats) %>%
    modify_column_hide(c(ci))
}


t5_total_anyHPV <- genotype %>%
  tbl_summary(
    by = any_hpv, 
    include = c(age_cat, M.STATUS, EDU.LEVEL, HIV, STI_RESULT2, 
                MSM.STATUS.x, S.WORKER, num_sex_partner, all_partner, 
                insert_partner, recep_partner, condom_partner), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "column"
  ) %>% 
  modify_column_hide(stat_1) %>% 
  modify_header(update = stat_2 ~ "**n(%), N = {n}**")


t5_cOR_anyHPV <- tbl_stack(list(
  glm_uni("any_hpv", "age_cat", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "M.STATUS", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "EDU.LEVEL", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "HIV", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "STI_RESULT2", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "MSM.STATUS.x", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "num_sex_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "all_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "insert_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "recep_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hpv", "condom_partner", genotype, "**cOR (95% CI)**")
))
t5_aOR_anyHPV <- glm_uni("any_hpv", 
                         paste0("age_cat + HIV + num_sex_partner"),
                         genotype, "**aOR (95% CI)**")
t5_anyHPV <- tbl_merge(list(t5_total_anyHPV, t5_cOR_anyHPV, t5_aOR_anyHPV), 
          tab_spanner = FALSE)



t5_total_anyHRhpv <- genotype %>%
  tbl_summary(
    by = any_hr_hpv, 
    include = c(age_cat, M.STATUS, EDU.LEVEL, HIV, STI_RESULT2, 
                MSM.STATUS.x, S.WORKER, num_sex_partner, all_partner, 
                insert_partner, recep_partner, condom_partner), 
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(
      all_continuous() ~ 1, 
      all_categorical() ~ c(0, 1)
    ),
    missing = "no", 
    percent = "column"
  ) %>% 
  modify_column_hide(stat_1) %>% 
  modify_header(update = stat_2 ~ "**n(%), N = {n}**")

t5_cOR_anyHRhpv <- tbl_stack(list(
  glm_uni("any_hr_hpv", "age_cat", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "M.STATUS", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "EDU.LEVEL", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "HIV", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "STI_RESULT2", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "MSM.STATUS.x", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "num_sex_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "all_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "insert_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "recep_partner", genotype, "**cOR (95% CI)**"),
  glm_uni("any_hr_hpv", "condom_partner", genotype, "**cOR (95% CI)**")
))
t5_aOR_anyHRhpv <- glm_uni("any_hr_hpv", 
                         paste0("age_cat + M.STATUS + EDU.LEVEL + HIV + STI_RESULT2 + ", 
                                "num_sex_partner + all_partner"),
                         genotype, "**aOR (95% CI)**")
t5_anyHRhpv <- tbl_merge(list(t5_total_anyHRhpv, t5_cOR_anyHRhpv, t5_aOR_anyHRhpv), 
                       tab_spanner = FALSE)

t5 <- tbl_merge(list(t5_anyHPV, t5_anyHRhpv), 
          tab_spanner = c("Any HPV", "Any HR-HPV")) %>%
  as_flex_table()  %>%
  set_caption("Regression analyses for HPV-associated risk factors")
t5 
        








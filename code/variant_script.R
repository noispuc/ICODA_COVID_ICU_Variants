#################################################################
###### Analysis of variants of concern in full covid cohort #####
#################################################################

# Libraries -----------------------------------------------------
library(tidyverse)
library(lubridate)
library(cowplot)
library(patchwork)
library(scales)
library(flextable)
library(gtsummary)
library(readxl)
library(performance)
library(parameters)
library(olsrr)
library(see)
library(lme4)
library(tidyverse)
library(flextable)
library(gtsummary)
library(parameters)
library(ggeffects)
library(emmeans)
library(rms)


# Data Input ----------------------------------------------------
df_covid_ext22 <- read.csv("df_icu_adm_MV_2022-05-31_COVID.csv")

colnames(df_covid_ext22)

### Descriptive analysis comparing dominant variants-------------

colnames(df_covid_ext22)

table0 <- df_covid_ext22 %>%
  filter((!is.na(mort60d))) %>%
  select(Age,
         idade_grupo,
         Gender,
         BMI,
         source_er,
         covid_primary,
         CharlsonComorbidityIndex,
         MFI_level,
         ChronicHealthStatusName,
         imunossupression,
         hypertension:obesity,
         dominant_variant) %>% 
  tbl_summary(by = dominant_variant,
              missing = "no",
              sort = all_categorical() ~ "frequency") %>%
  add_overall() %>%
  add_n() %>%
  add_p() %>%
  bold_labels() %>%
  as_flex_table()

table0

save_as_docx(table0, path = "table_baseline_3v_new_101022.docx")

table1 <- df_covid_ext22 %>%
  filter((!is.na(mort60d))) %>%
  select(Saps3Points,
         SofaScore,
         SOFARev,
         is_niv_24h,
         IsNeurologicalComaStuporObtundedDelirium,
         IsCardiovascularSepticShock,
         IsRespiratoryFailure:IsVasopressors,
         IsRenalReplacementTherapy,
         LowestGlasgowComaScale1h,
         ResourceIsMechanicalVentilation,
         ResourceIsRenalReplacementTherapy,
         ResourceIsVasopressors,
         ResourceIsTracheotomy,
         dominant_variant) %>% 
  tbl_summary(by = dominant_variant,
              missing = "no",
              sort = all_categorical() ~ "frequency") %>%
  add_overall() %>%
  add_n() %>%
  bold_labels() %>%
  add_p() %>%
  as_flex_table()

table1

save_as_docx(table1, path = "table_adm_comp_3v_new_101022.docx")

table2 <- df_covid_ext22 %>%
  filter((!is.na(mort60d))) %>%
  select(UnitDischargeCode,
         HospitalDischargeCode,
         UnitLengthStay,
         HospitalLengthStay,
         dominant_variant,
         first_MV_days,
         total_MV_days,
         los30d,
         los60d,
         mort30d,
         mort60d
  ) %>% 
  tbl_summary(by = dominant_variant,
              missing = "no",
              sort = all_categorical() ~ "frequency") %>%
  add_overall() %>%
  add_n() %>%
  bold_labels() %>%
  add_p() %>%
  bold_labels() %>%
  as_flex_table()

table2

save_as_docx(table2, path = "table_outcomes_3v_new_101022.docx")

# Estimating RCS knots for modelling -------------------------------------

knot_saps <- rcspline.eval(df_covid_ext22$Saps3Points, knots.only = TRUE)
knot_sofa <- rcspline.eval(df_covid_ext22$SofaScore, knots.only = TRUE)

#### MV model from "covid_ext_may22.R" -----------------------------------

mm_full <- glmer(mort60d ~ 1 + 
                   idade_grupo +
                   Gender +
                   dominant_variant +
                   covid_primary + 
                   source_er +
                   rcs(Saps3Points, parms = knot_saps) +
                   (1 | quarter) +
                   rcs(SofaScore, parms = knot_sofa) +
                   ChronicHealthStatusName +
                   MFI_level +
                   imunossupression +
                   hypertension +
                   cerebro_disease +
                   obesity +
                   IsRenalReplacementTherapy +
                   IsMechanicalVentilation +
                   IsVasopressors +
                   IsNeurologicalComaStuporObtundedDelirium +
                   (1 | UnitCode), 
                 data = df_covid_ext22, 
                 family = "binomial",
                 nAGQ = 0,
                 control = glmerControl(optCtrl = list(maxfun = 1e6)))


mm_full
summary(mm_full)
model_parameters(mm_full, exponentiate = TRUE,  details = TRUE)
check_model(mm_full)
model_performance(mm_full)

### Estimate marginal means for pairwise comparisons of variants -----------

vrt.emm <- emmeans(mm_full, params = c("knot_saps", "knot_sofa"),
                   specs = "dominant_variant",
                   rg.limit = 3000000,
                   type = "response",
                   adjust = "bonferroni")

pairs(vrt.emm, reverse = FALSE)
confint(vrt.emm)
confint(pairs(vrt.emm, reverse = FALSE))

### Estimate marginal means for mv model with interaction of ### -----------
### Variants and Mechanical Ventilation ### --------------------------------

mm_full_inter1 <- glmer(mort60d ~ 1 + 
                          idade_grupo +
                          Gender +
                          dominant_variant +
                          covid_primary + 
                          source_er +
                          rcs(Saps3Points, parms = knot_saps) +
                          (1 | quarter) +
                          rcs(SofaScore, parms = knot_sofa) +
                          ChronicHealthStatusName +
                          MFI_level +
                          imunossupression +
                          hypertension +
                          cerebro_disease +
                          obesity +
                          IsRenalReplacementTherapy +
                          IsMechanicalVentilation +
                          IsVasopressors +
                          IsNeurologicalComaStuporObtundedDelirium +
                          dominant_variant:IsMechanicalVentilation +
                          (1 | UnitCode), 
                        data = df_covid_ext22, 
                        family = "binomial",
                        nAGQ = 0,
                        control = glmerControl(optCtrl = list(maxfun = 1e6)))

model_parameters(mm_full_inter1, exponentiate = TRUE,  details = TRUE)

pp1 <-
  emmeans(
    mm_full_inter1,
    params = c("knot_saps", "knot_sofa"),
    pairwise ~ dominant_variant | IsMechanicalVentilation,
    rg.limit = 3000000,
    type = "response",
    adjust = "tukey",
    at = list(IsMechanicalVentilation = c(0, 1))
  )

summary(pp1)
confint(pp1$contrasts)
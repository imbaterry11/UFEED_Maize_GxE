#load libraries
require(chillR)
require(dormancyR)
require(pracma)
require(lubridate)
require(tidyr)
require(dplyr)
require(readr)
require(zoo)
require(fruclimadapt)
require(weathermetrics)
require(ggplot2)

trait_data = read_csv("1_Submission_Template_2024.csv")

##submission file prep

trait_data_submi = trait_data
y_pred_env = read_csv('yield_data_env_pred.csv')
y_pred_geno = read_csv('yield_data_geno_pred.csv')

trait_data_submi$Yield_Mg_ha = y_pred_env$Yield_Mg_ha
trait_data_submi$geno_impact = y_pred_geno$geno_impact

trait_data_submi= trait_data_submi %>%
  group_by(Env) %>%
  mutate(Yield_Mg_ha = mean(Yield_Mg_ha),
         Yield_Mg_ha = Yield_Mg_ha + geno_impact) %>%
  select(-geno_impact)

write_csv(trait_data_submi_3, file = 'trait_data_submi.csv')
# UFEED Maize GxE
This repo includes codes for generating the models needed for 2024 Maize_GxE competition. With the target of predicting 2024 maize yield data across hundreds of genotypes planted at various sites in U.S., this modeling approach is an example of combining Universal Feature Extraction from Environmental Data (UFEED) and Automated machine learning to predict 1) harvest date from each env (model 1); 2) mean yield per env (model 2); and 3) genotype effect and GxE effect on yield per env (model 3) **without any genotype data processing or genetic modeling**.
## All scripts and files needed for model generation and prediction are included in the main branch in this repo
### Follow this order to generate models and do model prediction
1) Run 'UFEED_GxE_train_data_processing.R' in R to generate input files for model training. <br />
   The files generated from this step are:
   * 'harvest_data_train.csv': input training file for the harvest date prediction model. Ready to use model can be found here:
   * 'yield_data_train_per_env.csv': input training file for the mean yield per env model and the genotype effect model. Ready to use model can be found here:
   * 'yield_data_test_per_env.csv': internal testing file for the mean yield per env model and the genotype effect model (internal testing was not included in the script). Ready to use model can be found here:
2) Run 'UFEED_models_training.ipynb' in Jupyter Notebook to generate the three models. <br />
   The folders generated from this step are:
   * 'maize_harvest_model': an AutoGluon model for the prediction of harvest date
   * 'maize_yield_model_env': an AutoGluon model for the prediction of mean yield per env
   * 'maize_yield_model_geno': an AutoGluon model for the prediction of genotype effect
3) Run 'UFEED_GxE_test_data_harvest_date_pred.R' in R to generate input files for harvest date prediction. <br />
   The file generated from this step is:
   * 'harvest_data_test.csv': input file for harvest date prediction
4) Run 'UFEED_models_harvest_date_prediction.ipynb' in Jupyter Notebook to predict harvest date. **There has to be a 'maize_harvest_model' folder containing the AutoGluon harvest date model in the same dir**<br />
   The files generated from this step are:
   * 'harvest_date_pred.csv'
   * 'pred_probs.csv'
5) Run 'UFEED_GxE_test_data_processing.R' in R to generate input files for yield prediction. <br />
   The file generated from this step is:
   * 'yield_data_predict.csv': input file for mean yield per env and genotype impact predictions
6) Run 'UFEED_models_yield_prediction.ipynb' in Jupyter Notebook to predict mean yield per env and genotype effect. **There has to be a 'maize_yield_model_env' folder containing the AutoGluon mean yield per env model and a 'maize_yield_model_geno' folder containing the AutoGluon genotype effect model in the same dir**<br />
   The files generated from this step are:
   * 'yield_data_env_pred.csv'
   * 'yield_data_geno_pred.csv'
8) Run 'UFEED_GxE_submission.R' in R to assemble files for final submission.<br />
   The file generated from this step is:
   * 'trait_data_submi.csv': ready-to-go file for submission

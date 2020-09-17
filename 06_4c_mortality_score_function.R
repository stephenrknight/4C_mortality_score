# ISARIC WHO CCP-UK study: 4C Mortality Score
# 4c mortality score function
# 06_4c_mortality_score_function.R
# Centre for Medical Informatics, Usher Institute, University of Edinburgh 2020

# 1. Apply variables changes across imputed datasets
# 2. glmnet using MICE datasets
# 3. Extract glmnet coefficients, combine and scale


# 4C mortality score function
isaric_lasso <- function(.data, age, sex, ethnicity, comorbid, rr, spo2, gcs, crp, bun,
                         output = c("vector", "components", "df_vector", "df_components"),
                         na_to_zeros = TRUE, all_na_rm = TRUE){
  .age = enquo(age)
  .sex = enquo(sex)
  .ethnicity = enquo(ethnicity)
  .comorbid = enquo(comorbid)
  .rr = enquo(rr)
  .spo2 = enquo(spo2)
  .gcs = enquo(gcs)
  .bun = enquo(bun)
  .crp = enquo(crp)
  
  out = .data %>%
    dplyr::mutate_at(vars(!! .age, !! .comorbid, !! .rr, !! .spo2, !! .gcs, !! .crp, !! .bun
    ), as.numeric) %>%
    dplyr::mutate_at(vars(!! .sex, !! .ethnicity), ~ as.character(.) %>% tolower() %>% trimws()) %>%
    
    mutate(
      isaric_lasso_age = case_when(
        !! .age < 50 ~ 0,
        !! .age < 60 ~ 2,
        !! .age < 70 ~ 4,
        !! .age < 80 ~ 6,
        !! .age >= 80 ~ 7,
        TRUE ~ NA_real_),
      
      isaric_lasso_sex = case_when(
        !! .sex == "male" ~ 1, # Changed from 0
        !! .sex == "female" ~ 0,
        TRUE ~ NA_real_),
      
      isaric_lasso_comorbid = case_when(
        !! .comorbid == 0 ~ 0,
        !! .comorbid == 1 ~ 1,
        !! .comorbid >= 2 ~ 2,
        TRUE ~ NA_real_),
      
      isaric_lasso_rr = case_when(
        !! .rr < 20 ~ 0,
        !! .rr < 30 ~ 1,
        !! .rr >= 30 ~ 2,
        TRUE ~ NA_real_),
      
      isaric_lasso_spo2 = case_when(
        !! .spo2 < 92 ~ 2,
        !! .spo2 >= 92 ~ 0,
        TRUE ~ NA_real_),
      
      isaric_lasso_gcs = case_when(
        !! .gcs <  15 ~ 2,
        !! .gcs == 15 ~ 0,
        TRUE ~ NA_real_),
      
      isaric_lasso_bun = case_when(
        !! .bun <= 7 ~ 0,
        !! .bun <= 14 ~ 2,
        !! .bun >  14 ~ 3,
        TRUE ~ NA_real_),
      
      isaric_lasso_crp = case_when(
        !! .crp < 50 ~ 0,
        !! .crp < 100 ~ 1,
        !! .crp >=100 ~ 2,
        TRUE ~ NA_real_)
    ) %>%
    mutate(
      isaric_lasso = rowSums(dplyr::select(., dplyr::starts_with("isaric_lasso_")),
                             na.rm = na_to_zeros)
    ) %>%
    
    {if(all_na_rm){
      dplyr::mutate(., isaric_lasso = dplyr::if_else(
        is.na(isaric_lasso_age) &
          is.na(isaric_lasso_sex) &
          is.na(isaric_lasso_comorbid) &
          is.na(isaric_lasso_rr) &
          is.na(isaric_lasso_spo2) &
          is.na(isaric_lasso_gcs) &
          is.na(isaric_lasso_bun) &
          is.na(isaric_lasso_crp), NA_real_, isaric_lasso))
    } else {
      .
    }}
  
  if(output == "vector"){
    out %>%
      pull(isaric_lasso)
  } else if(output == "components"){
    out %>%
      select(starts_with("isaric_lasso"))
  } else if(output == "df_vector"){
    out %>%
      pull(isaric_lasso) %>%
      bind_cols(.data, isaric_lasso = .)
  } else if(output == "df_components"){
    out
  }
}

# Prognostic index discrimination in mice derivation data ----------------------------------------
na_decision = FALSE   # Change to TRUE to impute 0 for NA in scoring models

sets_train %>% 
  complete("all") %>% 
  map(~ isaric_lasso(., age = age, sex = sex,  ethnicity = ethnicity_4levels, comorbid = no_comorbid, rr = rr_vsorres,
                     spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
                     bun = daily_bun_lborres, crp = daily_crp_lborres,
                     # plts = daily_plt_lborres, creat = daily_creat_lborres, nlr = NLR,
                     output = c("df_vector"), na_to_zeros = na_decision)
  ) %>% 
  map(~ roc(.$status, .$isaric_lasso) %>% 
        pROC::ci(method = "delong") # or bootstrap
  ) %>% 
  enframe() %>% 
  unnest_wider(value) %>% 
  summarise_all(mean)

# Prognostic index discrimination in mice validation data ----------------------------------------
sets_test %>% 
  complete("all") %>% 
  map(~ isaric_lasso(., age = age, sex = sex,  ethnicity = ethnicity_4levels, comorbid = no_comorbid, rr = rr_vsorres,
                     spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
                     bun = daily_bun_lborres, crp = daily_crp_lborres,
                     output = c("df_vector"), na_to_zeros = na_decision)
  ) %>% 
  map(~ roc(.$status, .$isaric_lasso) %>% 
        pROC::ci(method = "delong")
  ) %>% 
  enframe() %>% 
  unnest_wider(value) %>% 
  summarise_all(mean)





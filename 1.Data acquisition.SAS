*1-Extract covariates from the original data frame;

/* Define the library reference for the original data location */
LIBNAME ukb_covariates12 "U:\UK Biobank\OLINK\Full data set";

/* Extract specified covariates from the original dataset */
DATA mortality.ukb_covariates12 (KEEP = id age sex education smoking activity alcohol BMI hypertension diabetes dyslipidemia CVD cancer);
    SET ukb_covariates12.ukb675420; 
RUN;

*2-Extract metabolomic data from the original data frame;

/* Define the library reference for the outcomes data */
LIBNAME ukb_metabolites249 "U:\UK Biobank\OLINK\Full data set";

/* Extract metabolomics data from the original dataset */
DATA mortality.ukb_metabolites249 (KEEP = id Total_C HDL_C Total_TG Creatinine non_HDL_C Remnant_C VLDL_C Clinical_LDL_C LDL_C VLDL_TG LDL_TG HDL_TG Total_PL VLDL_PL LDL_PL HDL_PL Total_CE VLDL_CE LDL_CE HDL_CE Total_FC VLDL_FC LDL_FC HDL_FC Total_L VLDL_L LDL_L HDL_L Total_P VLDL_P LDL_P HDL_P VLDL_size LDL_size HDL_size Phosphoglyc TG_by_PG Cholines Phosphatidylc Sphingomyelins ApoB ApoA1 ApoB_by_ApoA1 Total_FA Unsaturation Omega_3 Omega_6 PUFA MUFA SFA LA DHA Omega_3_pct Omega_6_pct PUFA_pct MUFA_pct SFA_pct LA_pct DHA_pct PUFA_by_MUFA Omega_6_by_Omega_3 Ala Gln Gly His Total_BCAA Ile Leu Val Phe Tyr Glucose Lactate Pyruvate Citrate bOHbutyrate Acetate Acetoacetate Acetone Albumin GlycA XXL_VLDL_P XXL_VLDL_L XXL_VLDL_PL XXL_VLDL_C XXL_VLDL_CE XXL_VLDL_FC XXL_VLDL_TG XL_VLDL_P XL_VLDL_L XL_VLDL_PL XL_VLDL_C XL_VLDL_CE XL_VLDL_FC XL_VLDL_TG L_VLDL_P L_VLDL_L L_VLDL_PL L_VLDL_C L_VLDL_CE L_VLDL_FC L_VLDL_TG M_VLDL_P M_VLDL_L M_VLDL_PL M_VLDL_C M_VLDL_CE M_VLDL_FC M_VLDL_TG S_VLDL_P S_VLDL_L S_VLDL_PL S_VLDL_C S_VLDL_CE S_VLDL_FC S_VLDL_TG XS_VLDL_P XS_VLDL_L XS_VLDL_PL XS_VLDL_C XS_VLDL_CE XS_VLDL_FC XS_VLDL_TG IDL_P IDL_L IDL_PL IDL_C IDL_CE IDL_FC IDL_TG L_LDL_P L_LDL_L L_LDL_PL L_LDL_C L_LDL_CE L_LDL_FC L_LDL_TG M_LDL_P M_LDL_L M_LDL_PL M_LDL_C M_LDL_CE M_LDL_FC M_LDL_TG S_LDL_P S_LDL_L S_LDL_PL S_LDL_C S_LDL_CE S_LDL_FC S_LDL_TG XL_HDL_P XL_HDL_L XL_HDL_PL XL_HDL_C XL_HDL_CE XL_HDL_FC XL_HDL_TG L_HDL_P L_HDL_L L_HDL_PL L_HDL_C L_HDL_CE L_HDL_FC L_HDL_TG M_HDL_P M_HDL_L M_HDL_PL M_HDL_C M_HDL_CE M_HDL_FC M_HDL_TG S_HDL_P S_HDL_L S_HDL_PL S_HDL_C S_HDL_CE S_HDL_FC S_HDL_TG XXL_VLDL_PL_pct XXL_VLDL_C_pct XXL_VLDL_CE_pct XXL_VLDL_FC_pct XXL_VLDL_TG_pct XL_VLDL_PL_pct XL_VLDL_C_pct XL_VLDL_CE_pct XL_VLDL_FC_pct XL_VLDL_TG_pct L_VLDL_PL_pct L_VLDL_C_pct L_VLDL_CE_pct L_VLDL_FC_pct L_VLDL_TG_pct M_VLDL_PL_pct M_VLDL_C_pct M_VLDL_CE_pct M_VLDL_FC_pct M_VLDL_TG_pct S_VLDL_PL_pct S_VLDL_C_pct S_VLDL_CE_pct S_VLDL_FC_pct S_VLDL_TG_pct XS_VLDL_PL_pct XS_VLDL_C_pct XS_VLDL_CE_pct XS_VLDL_FC_pct XS_VLDL_TG_pct IDL_PL_pct IDL_C_pct IDL_CE_pct IDL_FC_pct IDL_TG_pct L_LDL_PL_pct L_LDL_C_pct L_LDL_CE_pct L_LDL_FC_pct L_LDL_TG_pct M_LDL_PL_pct M_LDL_C_pct M_LDL_CE_pct M_LDL_FC_pct M_LDL_TG_pct S_LDL_PL_pct S_LDL_C_pct S_LDL_CE_pct S_LDL_FC_pct S_LDL_TG_pct XL_HDL_PL_pct XL_HDL_C_pct XL_HDL_CE_pct XL_HDL_FC_pct XL_HDL_TG_pct L_HDL_PL_pct L_HDL_C_pct L_HDL_CE_pct L_HDL_FC_pct L_HDL_TG_pct M_HDL_PL_pct M_HDL_C_pct M_HDL_CE_pct M_HDL_FC_pct M_HDL_TG_pct S_HDL_PL_pct S_HDL_C_pct S_HDL_CE_pct S_HDL_FC_pct S_HDL_TG_pct);
    SET ukb_metabolites249.ukb675420;
RUN;

*3-Extract outcomes and follow-up time from the original data frame;

/* Define the library reference for the outcomes data */
LIBNAME ukb_outcomes "U:\UK Biobank\OLINK\Full data set";

/* Extract outcomes and follow-up time from the original dataset */
DATA mortality.ukb_outcomes (KEEP = id allcause_mortality followup_year);
    SET ukb_outcomes.ukb675420; /* Read the original dataset */
RUN;

*4-Merge above datasets;

PROC SORT DATA = mortality.ukb_covariates12; 
    BY id; 
RUN;

PROC SORT DATA = mortality.ukb_metabolites249; 
    BY id; 
RUN;

PROC SORT DATA = mortality.ukb_outcomes; 
    BY id; 
RUN;

/* Merge the datasets by id */
DATA mortality.ukb_fulldata;
    MERGE mortality.ukb_covariates12 
                  mortality.ukb_metabolites249 
                  mortality.ukb_outcomes;
    BY id;
RUN;

*5-Extract the dataset for analysis;

/* Process dataset mortality.ukb_fulldata */
DATA mortality.ukb_fulldata_final;
    SET mortality.ukb_fulldata;
    
    /* Retain rows where HDL Cholesterol is not missing */
    WHERE not missing(HDL_C); 
    
    /* Filter data for individuals aged between 50 and 69 */
    IF age >= 50 AND age <= 69;
    
    /* Create the variables for mortality after 10 years and 5 years */
    /* Mortality for all causes */
    IF followup_year > 10 THEN mortality10y = (allcause_mortality = 0);
    ELSE mortality10y = allcause_mortality;
    
    IF followup_year > 5 THEN mortality5y = (allcause_mortality = 0);
    ELSE mortality5y = allcause_mortality;
    
RUN;

*6-Repeat the above five steps to extract a comparable dataset from the ESTHER Study database, similar to the UK Biobank;

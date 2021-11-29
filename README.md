# Paradoxical
Code and data for: 
Paradoxical relationship between speed and accuracy in olfactory figure-background segregation.
Lior Lebovich, Michael Yunerman, Viviana Scaiewicz, Yonatan Loewenstein, and Dan Rokni


Raw data
Available in: all_mice_simplified_data_tables.mat
For data structure description see prepareRawData.m (line 3-14).


Script files:

prepareRawData.m
Coverts the relevant data in all_mice_simplified_data_tables.mat to goodAllMice_id_stim_rtNorm_choiceQ_0.csv
goodAllMice_id_stim_rtNorm_choiceQ_0.csv is used by all script below.

Figure2Figure3FigureS3.m
Plots Figure 2, Figure 3 and Suppplementary Figure 3.
(uses goodAllMice_id_stim_rtNorm_choiceQ_0.csv)

miceQuantile_regularModel.m
DDM fit of mice response data using Chi-square quantile optimization procedure and simulate data, based on DDM fit.
Outputs are saved to DDM\subj_params_miceq_chi.csv and DDM\simData_miceq_chi.csv, respectively.
(uses goodAllMice_id_stim_rtNorm_choiceQ_0.csv)

Figure4FigureS1FigureS2.m
Plots Figure 4,  Suppplementary Figure 1 and Suppplementary Figure 2.
(uses goodAllMice_id_stim_rtNorm_choiceQ_0.csv, DDM\subj_params_miceq_chi.csv and DDM\simData_miceq_chi.csv)



# Experiments & Real Data Applications



## Simulations

This folder contains code for implementing the simulation studies discussed in the paper. The correspondence between folders and figures/tables is summarized as follows:

```text
├── sim1-consistency
│   ├── DGPs
│   ├── binaryM		← Figure 3, 7
│   ├── continuousM	← Figure 4, 8
│   ├── multiM-d2	← Figure 5, 9
│   └── multiM-d4	← Figure 6, 10
├── sim2-weak-overlap	← Table 1, 8
│   ├── DGPs
│   ├── binaryM
│   ├── continuousM
│   └── multiM-d2
├── sim3-misspecification	← Table 2, 9
│   ├── DGPs
│   ├── binaryM
│   └── continuousM
├── sim4-crossfitting
│   ├── DGPs
│   ├── dense_forest	← Table 3, 10
│   └── sparse_forest	← Table 11, 12
├── sim5-sensitivity
│   ├── binaryY-saturated		|
│   ├── continuousY-binaryX		|	← Table 13
│   ├── continuousY-continuousX		|
│   ├── continuousY-continuousXM	← Table 4
│   ├── continuousY-continuousX_complex		← Table 14
│   
├── sim6-verma-efficiency
│   ├── DGPs
│   ├── R
│   ├── binaryM-binaryZ		← Figure 12
│   ├── binaryM-continuousZ		|
│   ├── binaryM-continuousZ-nomal01	| 	← Figure 13
│   ├── binaryM-continuousZ-nomal10	|
│   └── trash
├── sim7-nonlinearTMLE	← Table 6, 7
│   ├── DGPs
│   ├── binaryM		
│   ├── continuousM
│   └── multiM-d2
```



A summary of the functionality of commonly used files is provided below.

- joblist*.txt: This is the job file for simulation. Each line corresponding to one simulation. It is recommended to execute the job lists using parallel computing.

- write_job.R: This is the R script for producing the joblist*.txt files.

- main.R: Each line in the job list calls this main.R function to perform TMLE and one-step estimation. This file calls the 'fdtmle' package for estimation and save estimation results to the output folders, located under subfolders named after the estimators.

- organize_onestep.R: This file is used for organizing the output file from one-step estimators. It is called by the organize.txt file within each estimator folder.

- organize_onestep.R: This file is used for organizing the output file from TMLEs. It is called by the organize.txt file within each estimator folder.

- organize.txt: This file contains code for summarizing the files in the output folder. Run `bash organize.txt` in terminal to execute.

- plot.R: This is used for generating plots for sim1-consistency. This file calls plot-sub.R for generating smaller plots.

- plot-sub.R: This function is called by the plot.R for generating smaller plots.

- table.R: This is the R script used for generating tables in the paper.



## Real data application

```
B_PROUD
├── main_ordinal.R ← Produce estimates in Section 8
├── main_other.R ← Produce estimates in Appendix Section S8.1
├── organize_ordinal.R ← Organize estimates from main_ordinal.R
├── organize_other.R ← Organize estimates from main_other.R and create Appendix Table S12
├── write_job_ordinal.R
└── write_job_other.R
FSD
├── analysis_fsd.R ← Produce estimates in Appendix Section S8.2 
```



Utilizing our front-door estimation framework, we analyzed the impact of Mobile Stroke Unit (MSU) dispatch on patients’ functional outcomes using data from the [Berlin prehospital stroke care trial (B_PROUD)](ClinicalTrials.gov), under identifier: NCT02869386.



We also investigated how early academic achievement influences future annual income, using data from the [Life Course Study](https://services.fsd.tuni.fi/catalogue/FSD2076?tab=variables&lang=en&study_language=en).



Due to data sharing restrictions, we cannot provide direct access to the raw datasets used in this study. However, the dataset for the first application may be requested from the authors of the [original publication](https://journals.lww.com/epidem/abstract/2023/09000/the_effect_of_mobile_stroke_unit_care_on.14.aspx). The dataset for the second application is available through the [Finnish Social Science Data Archive](https://services.fsd.tuni.fi/catalogue/FSD2076?tab=variables&lang=en&study_language=en).




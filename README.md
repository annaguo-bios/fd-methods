# Experiments & Real Data Applications

This folder contains code for implementing the simulation studies and the real data application discussed in the paper.

## Simulations

The correspondence between folders and figures/tables is summarized as follows:

```text
├── sim1-consistency
│   ├── DGPs
│   ├── binaryM ← Figure S4, S8
│   ├── continuousM ← Figure S5, S9 
│   ├── multiM-d2 ← Figure S6, S10
│   └── multiM-d4 ← Figure S7, S11
├── sim2-weak-overlap ← Table 1, S4
│   ├── DGPs
│   ├── binaryM
│   ├── continuousM
│   └── multiM-d2
├── sim3-misspecification ← Table 2, S5
│   ├── DGPs
│   ├── binaryM
│   └── continuousM
├── sim4-crossfitting
│   ├── DGPs
│   ├── dense_forest ← Table S6, S7
│   └── sparse_forest ← Table S8, S9
├── sim5-sensitivity
│   ├── binaryY-saturated
│   ├── continuousY-binaryX ← Table S10
│   ├── continuousY-continuousX ← Table S10
│   ├── continuousY-continuousXM ← Table S10, S11
│   └── continuousY-continuousX_complex ← Table 3
├── sim6-verma-efficiency
│   ├── DGPs
│   ├── R
│   ├── binaryM-binaryZ-interAZ ← Figure S13
│   ├── binaryM-continuousZ-interAZ ← Figure S14
│   ├── binaryM-continuousZ-nonverma-estimator-interAZZ ← Figure S14
│   ├── binaryM-continuousZ-tn-interAZZ ← Figure S14
│   ├── binaryM-continuousZ-unif-interAZZ ← Figure S14
└── sim7-nonlinearTMLE ← Table S2, S3
    ├── DGPs
    ├── binaryM
    ├── continuousM
    └── multiM-d2
```



A summary of the functionality of commonly used files is provided below.

- joblist*.txt: This is the job file for simulation. Each line corresponding to one simulation. It is recommended to execute the job lists using parallel computing.

- write_job.R: This is the R script for producing the joblist*.txt files.

- main.R: Each line in the job list calls this main.R function to perform TMLE and one-step estimation. This file calls the 'fdcausal' package for estimation and save estimation results to the output folders, located under subfolders named after the estimators.

- organize_onestep.R: This file is used for organizing the output file from one-step estimators. It is called by the organize.txt file within each estimator folder.

- organize_TMLE.R: This file is used for organizing the output file from TMLEs. It is called by the organize.txt file within each estimator folder.

- organize.txt: This file contains code for summarizing the files in the output folder. Run `bash organize.txt` in terminal to execute.

- plot.R: This is used for generating plots for sim1-consistency. This file calls plot-sub.R for generating smaller plots.

- plot-sub.R: This function is called by the plot.R for generating sub plots.

- table.R: This is the R script used for generating tables in the paper.



## Real data application

The correspondence between folders and tables is summarized as follows:

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

We present two real-data applications to illustrate the practical use of our proposed front-door estimation framework:

- The first application uses data from the [Berlin Prehospital Stroke Care Trial (B_PROUD)](https://clinicaltrials.gov/study/NCT02869386) (ClinicalTrials.gov identifier: NCT02869386).
- The second application uses data from the [Life Course Study](https://services.fsd.tuni.fi/catalogue/FSD2076?tab=variables&lang=en&study_language=en).

### Data Availability

Due to data sharing restrictions, we cannot provide direct access to the original datasets used in this study.

- For the B_PROUD application, the dataset may be requested from the authors of the [original publication](https://journals.lww.com/epidem/abstract/2023/09000/the_effect_of_mobile_stroke_unit_care_on.14.aspx).
- For the Life Course Study, the dataset is available through the [Finnish Social Science Data Archive](https://services.fsd.tuni.fi/catalogue/FSD2076?tab=variables&lang=en&study_language=en).

### Synthetic Data

To facilitate reproducibility, we provide synthetic datasets for both applications:

- Each dataset is named `synthetic_data.csv`.
- The files are located in the `B_PROUD/` and `FSD/` directories, respectively.

These synthetic datasets are generated to preserve the structure of the original data and can be used to run and test the code in this repository.
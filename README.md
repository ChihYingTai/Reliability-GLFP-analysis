# Reliability-GLFP-analysis
Reliability inference in GLFP models based on EM algorithm with related application

Data Folder: 

There are 2 files in this folder. 

1. Gateoxide.cdv: exact failuretime information with columns:
  - endtime: the observed failure time or censored time of the product. 
  - censored: 1 for right-censored observations, 0 for failures.

2. CB_data.csv: interval censored failure time information with columns: 
  - starttime: lower bound of the time interval. 
  - endtime: upper bound of the time interval. 
  - count: number of failures occurring within the interval. 
  - censored: 1 for right-censored observations, 0 for failures.


Exact_failure Folder:

There are 11 files in this folder. 

1. Data_func: Function for generating exact failure time data.

2. EM_opt: EM algorithm implementation.

3. MLE_opt: MLE implementation.

4. FIM_complete_ob: Compute the observed FIM using the complete likelihood.

5. FIM_missing_ob: Compute the observed FIM using missing information.

6. FIM_SE_CI_tp_func: Estimate quantiles and corresponding SEs for failure modes 1, 2, and the system.

7. Predict_class: Predict failure mode and defectiveness; includes confusion matrix if true labels are provided.

8. Probability_plot: Generate probability plots for both exact failure data and interval-censored data.

9. Hazard_plot: Generate hazard plots for both exact failure data and interval-censored data.


Interval_censored Folder:

There are 13 files in this folder. 

1. EM_interval_opt: EM algorithm implementation for interval-censored data.

2. MLE_interval_opt: MLE implementation for interval-censored data.

3. FIM_interval_complete_ob:  Computes the observed FIM using the complete likelihood for interval-censored data.

4. FIM_interval_missing_ob: Computes the observed FIM using the missing information for interval-censored data.

5. Predict_interval_class: Predict failure mode and defectiveness for interval-censored data; includes confusion matrix if true labels are provided.









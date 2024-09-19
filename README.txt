This folder consists of the R codes for the method used to perform the simulation studies and real data analysis in the manuscript. In this README file, I will explain each R file.

###########################################################

1. FaFQL_est.R: the R code used to calculate the estimation for the parameters of the functional quantile linear regression with hidden factors (FaFQL). The initial output is used as the starting point of the iteration algorithm for computing the estimations of the proposed method.

2. FaFQL_fact.R: the R code used to identify the optimal factor for the proposed model when the true factor number is not given. 

3. FQLM_est.R: the R code used to calculate the conventional functional quantile linear model (FQLM) without considering the hidden factors in the proposed model.

4. FaFQL_eval.R: the R code used to evaluate the the performance of estimates in proposed FaFQL model and conventional functional linear quantile model (FQLM) without considering the hidden factors.

5. Sam_data.R: the R code used to generate data for the simulation studies presented in the manuscript.

6. Real_data.R: the R code used to conduct real data analysis presented in the manuscript.
###########################################################
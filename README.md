# High-Dimensional Model Averaging via Cross-validation
To reproduce the data analysis results in the paper, simply run the files: Simulation_Lin.R, Simulation_Log.R, Simulation_Inference_Lin.R, Simulation_Inference_Log.R and Real-Data.R.

Here are some descriptions for other R files.
(1) AnL: Leave-one-out cross-validation model averaging proposed by Ando and Li.
(2) PLASSO: Post Lasso (Refit an OLS linear model (or logistic regression model with Lasso penalty) using the variables selected by Lasso).
(3) PMA: Parsimonious model averaging in linear regression.
(4) HDMA: Our proposed model averaging approach. (Two loss functions: square loss in linear regression and cross-entropy loss in logistic regression; Three penalties: Lasso, SCAD and MCP; Two methods for solving simplex-constrained optimization problem: FGMA and GMA)
(5) Bootstrap: Implementing Gaussian multiplier bootstrap to construct simultaneous confidence intervals.
(6) Fun: Some functions we used for generating data and simple calculations.

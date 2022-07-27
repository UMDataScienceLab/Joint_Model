# Joint_Model
Joint Models for Event Prediction from Time Series and Survival Data, Technometrics 2020.

## Introduction and Contributions

* We present a nonparametric prognostic framework for individualized event prediction based on joint modeling of both time series and time-to-event data. 
* Our approach exploits a multivariate Gaussian convolution process (MGCP) to model the evolution of time series signals and a Cox model to map time-to-event data with time series data modeled through the MGCP. Taking advantage of the unique structure imposed by convolved processes, we provide a variational inference framework to simultaneously estimate parameters in the joint MGCP-Cox model. This significantly reduces computational complexity and safeguards against model overfitting. 
* Experiments on synthetic and real world data show that the proposed framework outperforms state-of-the art approaches built on two-stage inference and strong parametric assumptions. 

## Prerequisite

* [R](https://www.r-project.org/)

## Code:

* joint_estimation.R:

This is the main file that will learn a multivariate Gaussian process based joint model.
 
* new_rejectSam.R and R_functions_new_new.r

Those files contain auxillary functions such as rejection sampling, simulation data generation process, etc.

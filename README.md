# BAMVALS
Bayesian Framework for a MIMO Volterra Tensor Network

## Abstract
This paper proposes a Bayesian Volterra tensor network (TN) to solve high-order discrete nonlinear multiple-input multiple-output (MIMO) Volterra system identification problems. Using a low-rank tensor network to compress all Volterra kernels at once, we avoid the exponential growth of monomials with respect to the order of the Volterra kernel. Our contribution is to introduce a Bayesian framework for the low-rank Volterra TN. Compared to the least squares solution for Volterra TNs, we include prior assumptions explicitly in the model. In particular, we show for the first time how a zero-mean prior with diagonal covariance matrix corresponds to implementing a Tikhonov regularization for the MIMO Volterra TN. Furthermore, adopting a Bayesian viewpoint enables simulations with Bayesian uncertainty bounds based on noise and prior assumptions. In addition, we demonstrate via numerical experiments how Tikhonov regularization prevents overfitting in the case of higher-rank TNs.




## Before running the code 

To run the code, execute the script BAMVALS.m
For the code to run you must download the following files:

https://github.com/kbatseli/MVMALS [1]
https://data.4tu.nl/articles/_/12960104 [2]


[1] Batselier, Kim, Zhongming Chen, and Ngai Wong. "Tensor Network alternating linear scheme for MIMO Volterra system identification." Automatica 84 (2017): 26-35.
[2] Schoukens, Maarten, et al. "Cascaded tanks benchmark combining soft and hard nonlinearities." Workshop on nonlinear system identification benchmarks. 2016.

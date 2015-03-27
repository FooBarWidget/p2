This is an implementation of the P^2 algorithm as described in the paper "The P^2 Algorithm for Dynamic Calculation of Quantiles and Histograms Without Storing Observations" by Raj Jain and Imrich Chlamtac. It's based on [Aaron Small's implementation](https://github.com/absmall/p2).

I wrote this as part of my research into algorithms for computing quantiles without storing every observations. Other interesting resources:

 * [Algorithm to dynamically monitor quantiles](https://stats.stackexchange.com/questions/7959/algorithm-to-dynamically-monitor-quantiles/14937) on StackOverflow
 * "Monitoring Networked Applications With Incremental Quantile Estimation" by Chambers et al describes an alternative algorithm, one which is specifically designed for distributed data collection and aggregation. An implementation of this algorithm ("IQAgent") is described in "Numerical Recipes" by Press et al, section 8.5.2.
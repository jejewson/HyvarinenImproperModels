# HyvarinenImproperModels

HyvarinenImproperModels contains the *R* scripts and *stan* software to reproduce the experiments from the paper "General Bayesian Loss Function Selection and the use of Improper Models" J. Jewson and D. Rossell (2021). The article and supplementary material can be found at *https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12553*.

The repository contains the following .Rmd files:
```
Marginal_kappa.Rmd, 
Gaussian_Experiments.Rmd, 
TGFB_data.RMD, 
DLD_data.Rmd, 
KernelDensityEstimation.Rmd,
TsallisScore_ImproperModel.Rmd,
SeparatedMultiModal_experiments.Rmd
```
which respectively reproduces the experiments from Section 5.1, 5.2, 5.3.1, 5.3.2, 6.1, A.3 and A.5.3  of the paper. 

The *data* folder contain the TGF-$\beta$ and DLD datasets. The *R* folder contains the functions required to analytically evaluate the Hessian of the Gaussian model and Tukey's loss (under the smooth approximation to the absolute loss), as well as specify the non-local priors according to Section A.3.4 and evaluate kernel density estimates. The *stan* folder contains the necessary .stan files to calculate MAP estimates (and sample) form the $\mathcal{H}$ posteriors for the Gaussian model, Tukey's loss and the kernel density estimate.

### Contact Information and Acknowledgments

HyvarinenImproperModels was developed by Jack Jewson, Universitat Pompeu Fabra, Barcelona and The Barcelona School of Economics (jack.jewson@upf.edu). 

The project was partially funded by the Ayudas Fundación BBVA a Equipos de Investigación Cientifica 2017 and Government of Spain's Plan Nacional PGC2018-101643-B-I00 grants. 






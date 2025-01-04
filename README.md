## Predicting fine-scale spatiotemporal persistence of rare species: an example with North Atlantic right whales   

#### Jamie L. Brusa, Daniel W. Linden, Meghan P. Gahm, Eric M. Patterson, Caroline P. Good, Daniel E. Pendleton, Jason J. Roberts, Timothy V. N. Cole, Danielle M. Cholewiak  

##### Please contact the first author for questions about the code or data: Jamie Brusa (jlbwcc@gmail.com)

DOI: Note, this work is in the preparation phase an not yet submitted for publication (and this repo is currently set to private)
_______________________________________________________________________________________

## Abstract

Knowledge of where, when, and for how long rare species occur is essential for effective conservation and management. The use of species occurrence data to predict persistence (in space and time) may be possible with extensive survey effort and complex modeling, but such requirements present challenges. Here, we present a simple approach for estimating wildlife spatiotemporal persistence using a data filter and bootstrap resampling. Our method avoids the need for mechanistic models of animal behavior and habitat selection while accommodating survey effort that may, at times, be sparse. We illustrate the approach using spatiotemporal data from vessel-based and aerial surveys of the critically endangered North Atlantic right whale (Eubalaena glacialis), a rare and cryptic species that endures various anthropogenic disturbances and motivates controversial management actions. Further, persistence of this species has material effects on a diverse group of people, including commercial and recreational ocean users. Our analyses suggested that persistence probabilities of North Atlantic right whales varied across time and space with patterns that could be useful for guiding management measures. Our approach to estimating spatiotemporal persistence can be applied to any species with sufficient survey data and could facilitate dynamic management practices that more effectively target conservation efforts in time and space. The proposed method is especially helpful for rapid decision making when more complex models might require more time to update.


### Required Packages and Versions Used
####For data filtering and bootstrap resampling
tidyverse
here
sf
sp
stringer
####Additional package for testing for out-of-sample predictive ability
caret
####Additional packages for modelling effort and detections
raster
rjags
runjags
coda
ggmcmc
mcmcVIS
arm
loo


### How to Use this Repository
The data filtering and analysis code can be used to estimate the probability of detecting an individual of a particular species given a recent previous detection of the same species in the same geographic location for visual detection data of wildlife species. We provide code for data filtering and bootstrap resampling and for simulations (to test for out-of-sample predictive ability). Additionally, we also provide code described in Appendix 2 and Appendix 3 of the accompanying publication for modeling the relationship between survey effort and redetections as well as for incorporating prior information using a similar approach to our resampling protocol but using Bayes' Theorem. The code is versitile to be adapted to data for other studies that might benefit from our presented protocol.

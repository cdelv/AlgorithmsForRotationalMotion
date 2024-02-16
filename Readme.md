# SPIRAL: An Efficient Algorithm for the Integration of the Equation of Rotational Motion

This repo contains the software used to compare multiple rotation algorithms and the implementation of SPIRAL, a new algorithm that vastly outperforms traditional ones.

## Abstract

We introduce SPIRAL, a third-order integration algorithm for the rotational motion of extended bodies. It requires only one force calculation per time step, does not require quaternion normalization at each time step, and can be formulated for both leapfrog and synchronous integration schemes, making it compatible with many particle simulation codes. The stability and precision of SPIRAL exceed those of state-of-the-art algorithms currently used in popular DEM codes such as Yade, MercuryDPM, LIGGGHTS, PFC, and more, at only slightly higher computational cost. Also, beyond DEM, we see potential applications in all numerical simulations that involve the 3D rotation of extended bodies.

## Citation

```
@article{SPIRAL,
title = {SPIRAL: An efficient algorithm for the integration of the equation of rotational motion},
journal = {Computer Physics Communications},
volume = {297},
pages = {109077},
year = {2024},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2023.109077},
url = {https://www.sciencedirect.com/science/article/pii/S0010465523004228},
author = Carlos Andr\'es {del Valle} and Vasileios Angelidakis and Sudeshna Roy and Jos\'e Daniel Mu\~noz and Thorsten P\"oschel},
}
```

You can find the pre-print in Arxiv. 

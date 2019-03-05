# PRS_with_DCFTP
This repository consists of python codes for perfect sampling using partial rejection sampling (PRS) proposed by Moka and Kroese (2019). This method combines discrete PRS of Guo et al. (2017) and dominated coupling from the past (DCFTP). Details of the python codes are given below:


  1. PRS_PSM.py: generates samples from penetrable spheres mixture model on [0,1]^2, which is absolutely continuous with a homogeneous Poisson point process with intensity 4*(kappa_1 + kappa_2)/(pi*IntRange*IntRange). Where 4*kappa_i/(pi*IntRange) is the intensity of type_i points, and IntRange is the interaction range. On each cell, perfect samples generated using dominated CFTP.

  2. PRS_Strauss.py: generates samples from Strauss process on [0,1]^2, which is absolutely continuous with a homogeneous Poisson point process with intensity 4*kappa_0/(pi*IntRange*IntRange), where IntRange is the interaction range. On each cell, perfect samples generated using Huber's dominated CFTP.

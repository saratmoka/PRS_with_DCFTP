# PRS_with_DCFTP
PRS_PSM.py and PRS_Strauss.py are based on the partial rejection sampling (PRS) by Moka and Kroese (2019). PRS_GJ_HardCore.py is based on the PRS by Guo and Jerrum (2018). More description of these python codes are given below:

  1. PRS_PSM.py: generates samples from penetrable spheres mixture model on [0,1]^2, which is absolutely continuous with a homogeneous Poisson point process with intensity 4*(kappa_1 + kappa_2)/(pi*IntRange*IntRange). Where 4*kappa_i/(pi*IntRange) is the intensity of type_i points, and IntRange is the interaction range. On each cell, perfect samples generated using dominated CFTP.

  2. PRS_Strauss.py: generates samples from Strauss process on [0,1]^2 with parameter gamma, which is absolutely continuous with a homogeneous Poisson point process with intensity 4*kappa_0/(pi*IntRange*IntRange), where IntRange is the interaction range. On each cell, perfect samples generated using Huber's dominated CFTP. Strauss process becomes hard-core process when gamma = 0.
  
  3. PRS_GJ_HardCore.py: generates samples from Hard_core process on [0,1]^2 with interaction range IntRange and absolutely continuous with respect to a Poisson point process with intensity 4*kappa_0/(pi*IntRange*IntRange).

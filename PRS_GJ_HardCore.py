# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:55:29 2018

@author: uqsbabum
"""
def reset():
    try:
        get_ipython().magic('reset -sf')  #analysis:ignore
    except NameError:
        pass
reset()

# Program to generate hard-core process on [0,1]^2 with Poisson intensity beta and hard-core distance HCdist.
# It is based on the partial rejection sampling by Guo and Jerrum, 2018.


import numpy as np
import time
#import matplotlib.pyplot as plt
import datetime

## Euclidean distance between two points on R^2 
def euclidean2dim(yy,xx):
    return np.sqrt((xx[0] - yy[0])**2 + (xx[1] - yy[1])**2)


## Subroutine for generating n points on a circle with radius rad and ceter (x,y)
def PointsOnCircle(n, rad, x, y):
    Points = []
    for _ in range(n):
        r = rad*np.sqrt(np.random.random_sample())
        theta = 2*np.pi*np.random.random_sample()
        Points.append((x + r*np.cos(theta), y + r*np.cos(theta))) 
    return Points    


Beta_0 = 0.1


Est_time = 0
Exp_while_loops = 0
npoints = 0
Exp_total_pts = 0
HCdist = 1/100   #HCdist/2 is the radius of each sphere.
Beta = np.divide(4*Beta_0, np.pi*HCdist*HCdist)
Beta_circ = Beta*np.pi*HCdist*HCdist


Itot = 1000 ## No of perfect samples generated


print('\n Intensity = ', Beta, '\n Beta_0:', Beta_0, '\n Hard-core distance = ', HCdist,'\n # iterations = ', Itot)
print('Program Starting time: ', datetime.datetime.now().time())
print('++++++++++ GJ Alg for Hard Core Process +++++++++++++++') 

np.random.seed(0)

for b in range(Itot):
    print('---------- GJ PRS Iteration: ', b+1, '---------')
    N = np.random.poisson(Beta)
    Disks = [(np.random.random_sample(), np.random.random_sample()) for _ in range(N)]

    count_while = 0
    while True:
        count_while = count_while + 1
        
        N = len(Disks)
        Bad_ind = []
        
        stime = time.time()
        for n in range(N):
            Bad = False
            for m in range(N):
                if m != n and euclidean2dim(Disks[n], Disks[m]) < HCdist:
                    Bad = True
                    break
    
            if Bad == True:
                Bad_ind.append(n)
        print('Time 1:', time.time() - stime)
        len_bad = len(Bad_ind)
        
        if len_bad == 0:
            break
        else:
            Ps = PointsOnCircle(np.random.poisson(Beta_circ), HCdist, Disks[Bad_ind[0]][0], Disks[Bad_ind[0]][1])
            for i in range(len_bad - 1):
                PP_circ = PointsOnCircle(np.random.poisson(Beta_circ), HCdist, Disks[Bad_ind[i+1]][0], Disks[Bad_ind[i+1]][1])
                for x in PP_circ:
                    Inside = False
                    for j in range(i+1):
                        if euclidean2dim(x, Disks[Bad_ind[j]]) < HCdist:
                            Inside = True
                            break
                    if Inside == False:
                        Ps.append(x)
                        
            Disks_good = [Disks[n] for n in range(N) if (n not in Bad_ind)]
            Disks = Disks_good + Ps
 
    npoints = np.divide(b, b + 1)*npoints + np.divide(len(Disks), b + 1)
    print('exp no of points :', npoints)
    Exp_while_loops = np.divide(b, b + 1)*Exp_while_loops + np.divide(count_while, b + 1)
    print('Expected while loops :', Exp_while_loops)

print('++++++++++ GJ Alg for Hard Core Process +++++++++++++++') 
print('Parmeters: \n Intensity = ', Beta, '\n Beta_0:', Beta_0, '\n Hard-core distance = ', HCdist,'\n # iterations = ', Itot)  
print('\nResults:')
print('  Estimated # points:', npoints)
print('  Expected # while loops:', Exp_while_loops)
print('Program Ending time: ', datetime.datetime.now().time())
       




     
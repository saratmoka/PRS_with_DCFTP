# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:55:29 2018

@author: uqsbabum
"""

# Program to generate Strauss process on [0,1]^2 with parameter Gamma and abosolutely continuous with respect to 
# the Poisson point process with intensity Kappa = Kappa_0/(pi*IntRange*IntRange), where IntRange is 
# the interaction range of the model.


def reset():
    try:
        get_ipython().magic('reset -sf')  #analysis:ignore
    except NameError:
        pass
reset()


import numpy as np
import time
import datetime
import math




###### Euclidean distance between two points on R^2 
def euclidean2dim(yy, xx):
    return np.sqrt((xx[0] - yy[0])**2 + (xx[1] - yy[1])**2)

###### Returns True if there is a point in A within IntRange froma point in B. Otherwise, returns False.
def Overlap(A, B, IntRange):
    for x in A:
        for y in B:
            if euclidean2dim(x, y) < IntRange:
                return True
    return False  

###### Subroutine to count the number of overlaps the points in A have with the points in B
def OverlapCount(A, B, IntRange):
    count = 0
    for x in A:
        for y in B:
            if  euclidean2dim(x, y) < IntRange:
                count = count + 1
            
    return count

###### Subroutine to find the indecies of the circles in A that are blocking the new circle NC. Used in dCFTP function. See Huber's swap based dCFTP for details
def BlockInd(A, NC, IntRange):
    Ind = []
    
    for i in range(len(A)):
        if euclidean2dim(A[i][0:2], NC[0:2]) < IntRange and A[i][2] == 0:
            Ind.append(i)

    return Ind


###### Huber's dominated CFTP method for generating a perfect sample of Strauss process on Cell numbered 'Cell_no' with parameter gamma and intensity Kappa_cell
def dCFTP_swap(Cell_no, Kappa, gamma, IntRange):
    x_cord = Cell_no%K
    y_cord = int(Cell_no/K)
    Cell_x_len = min(IntRange, 1 - x_cord*IntRange)
    Cell_y_len = min(IntRange, 1 - y_cord*IntRange)

    Kappa_cell = Kappa*Cell_x_len*Cell_y_len
    
    D = [] 
    events = [] #Records the events happened
    D.append([(Cell_x_len*np.random.random_sample(), Cell_y_len*np.random.random_sample(), np.random.binomial(1, gamma)) for i in range(np.random.poisson(Kappa_cell))]) # Intial state 

    D.append(D[0][:])
    no_cirD = len(D[0]) # Number of disks at time 0

    if  no_cirD == 0 or np.random.binomial(1, np.divide(Kappa_cell, Kappa_cell + no_cirD)) == 1: #True if it is an arrival

        NewCen = (Cell_x_len*np.random.random_sample(), Cell_y_len*np.random.random_sample(), np.random.binomial(1, gamma))

        D[1].append(NewCen)
        no_cirD = no_cirD + 1
        events.append((-1, NewCen))

        
    else:
        k = np.random.randint(no_cirD) ## Index of the departuring element
        vv = D[1].pop(k)
        events.append((1,vv))
        no_cirD = no_cirD - 1 

    n = 1

    while True: ## This while Loop ends only when there is a coalscence
        L = [] # initialization of Lower bounding process
        U = D[n][:] # initialization of Upper bounding process

        for j in range(n):
            if events[n-j-1][0] == -1: ## It is a deprture in the forward chain (since it is an arrival backward in D)
                try:
                    v_index = L.index(events[n - j - 1][1])
                except:
                    v_index = -1
                    
                if v_index != -1:
                    del L[v_index]

                
                try:
                    v_index = U.index(events[n - j - 1][1])
                except:
                    v_index = -1
                
                if v_index != -1:
                    del U[v_index]
                
            else:
                Block_Ind_U = BlockInd(U, events[n - j - 1][1], IntRange)
                Block_Ind_L = BlockInd(L, events[n - j - 1][1], IntRange)
                
                if len(U) == 0 or len(Block_Ind_U) == 0:
                    L.append(events[n - j - 1][1])
                    U.append(events[n - j - 1][1])
                    
                elif len(Block_Ind_U) == 1:

                    try:
                        v_index = L.index(U[Block_Ind_U[0]])
                    except:
                        v_index = -1
                        
                    if v_index != -1:
                        del L[v_index]                        

                    del U[Block_Ind_U[0]]                        
                    
                    L.append(events[n - j - 1][1])
                    U.append(events[n - j - 1][1])                        
                    
                elif len(Block_Ind_L) == 0:
                    U.append(events[n - j - 1][1])
                    
                    
                elif len(Block_Ind_L) == 1:

                    del L[Block_Ind_L[0]]                        
                    U.append(events[n - j - 1][1])

            if (set(L) <= set(U)) == False:
                print("L is not a subset of U")
                exit()
        
        if len(L) == len(U):
            break

        for j in range(n): #Loop for doubling the length of the dominated process
            D.append(D[n+j][:])
            
            if  no_cirD == 0 or np.random.binomial(1, np.divide(Kappa_cell, Kappa_cell + no_cirD)) == 1: # True if it is an arrival
                NewCen = (Cell_x_len*np.random.random_sample(), Cell_y_len*np.random.random_sample(), np.random.binomial(1, gamma))
                D[n+j+1].append(NewCen)
                no_cirD = no_cirD + 1
                events.append((-1, NewCen))
            else:
                k = np.random.randint(0, no_cirD) ## Index of the departuring element
                vv = D[n+j+1].pop(k)
                events.append((1,vv))                    
                no_cirD = no_cirD - 1 
        
        n = 2*n            

    State = []
    for pt in L:
        State.append((pt[0] + x_cord*IntRange, pt[1] + y_cord*IntRange))
        
    return State


###### Node of the dependency graph
class GNode:
    def __init__(self, initdata):
        self.event = initdata
        self.uniform = 0
        self.next = None
        
###### Node for keeping a neighbor of a given node of the dependency graph
class NNode:
    def __init__(self, initdata):
        self.neighbor = initdata
        self.next = None

##### Generating the nodes of the dependency graph    
def generateVertices(K):
    stime = time.time()
    head = None
    Vert = set([])
    for x in range(K*K):
        x0 = x%K
        x1 = int(x/K)
        for y0 in range(max(0, x0-1), min(K, x0+2)):
            for y1 in range(max(0, x1-1), min(K, x1+2)):        
                Vert.add(frozenset((x, y1*K+y0)))
                
        Vert.remove(frozenset((x,x)))
        
    for a in Vert:
        pt = GNode(list(a))
        pt.next = head
        head = pt
        
    print("\n Vertices of the dependency graph is generated. Time taken: ", int(time.time() - stime), "secs")
                
    return head


###### Identifying the neighbors of each event in the dependency graph
def generateNeighbors(K):
    Neighbors = [None for _ in range(K*K)]
    for x in range(K*K):
        x0 = x%K
        x1 = int(x/K)
        Nb = [j*K+i for i in range(max(0, x0-1), min(K, x0+2)) for j in range(max(0, x1-1), min(K, x1+2))]

        Nb.remove(x)
        
        for y in Nb:
            head = NNode(y)
            head.next = Neighbors[x]
            Neighbors[x] = head
    return Neighbors


## New Partial Rejection Sampling Alg
def straussNPRS(K, Kappa, gamma, IntRange, Graph, Neighbors):
    
    Disks = [[] for _ in range(K*K)] 
    BadCells = [False for _ in range(K*K)]   
    
    ## Intial state generation
    for x in range(K*K):  # x is the cell number
        Disks[x] = dCFTP_swap(x, Kappa, gamma, IntRange)
        
    ## Generate uniform random variables associated with each event/edge 
    pt = Graph
    while pt != None:
        pt.uniform = np.random.random_sample()
        pt = pt.next
    
    # Start of the while loop of Partial Rejection Sampling Method
    count_while = 0
    while True:
        count_while = count_while + 1
    
        ## Identifying the bad events
        stop = True
        pt = Graph

        while pt != None:
            (x,y) = list(pt.event)

            count_overlaps = OverlapCount(Disks[x], Disks[y], IntRange)
                
            if pt.uniform > np.power(gamma, count_overlaps):
                stop = False
                BadCells[x] = True
                BadCells[y] = True

            pt = pt.next
    
        ## PRS is terminated if there is no bad event
        if stop:
            break

##### Construction of resampling cells
        ResCells = BadCells.copy()
        stop = False

        while not stop:
            stop = True
            for x in range(K*K):
                if ResCells[x]:
                    current = Neighbors[x]
                    while current != None:
                        y =  current.neighbor
                        if (not ResCells[y]) and (len(Disks[y]) > 0):
                            stop = False
                            ResCells[y] = True
                        current = current.next
                        
        pt = Graph
        while pt != None:
            (x,y) = list(pt.event)
            if ResCells[x] == True or ResCells[y] == True:
                pt.uniform = np.random.random_sample()
            pt = pt.next                        
                        
        for x in range(K*K):
            if len(Disks[x]) and ResCells[x]:
                current = Neighbors[x]
                while current != None:
                    ResCells[current.neighbor] = True
                    current = current.next
                        
        
        for x in range(K*K):
            if ResCells[x]:
                Disks[x] = dCFTP_swap(x, Kappa, gamma, IntRange)                
                        
            BadCells[x] = False
                

    ## Conclusions after the iteration                       
    npoints = 0
    for x in range(K*K):
        npoints = npoints + len(Disks[x])
    result = (npoints, count_while)
    
    return result

###### Model parameters
Kappa_0 = 0.1
gamma = 0.0
IntRange = 1/100

###### Partition parameters
K = math.ceil(1/IntRange)
Kappa = np.divide(4*Kappa_0, np.pi*IntRange*IntRange)


Itot = 1000 ### Number of samples generated

###### Estimation parameters
Est_time = 0
Exp_while_loops = 0
Exp_total_pts = 0
npoints = 0



print('Kappa_0 = ', Kappa_0, ', Gamma = ', gamma, '\n Interaction Range = ', IntRange,'\n # iterations = ', Itot)
print('No of cells = ', K*K)

np.random.seed(0)
print('Program Starting time: ', datetime.datetime.now().time())
print('++++++++++ New PRS Alg for Strauss Process +++++++++++++++')




Graph = generateVertices(K)
Neighbors = generateNeighbors(K)

for b in range(Itot):
    print('---------- New PRS Iteration: ', b+1, '---------')
    
    
    (count_pts, count_while) = straussNPRS(K, Kappa, gamma, IntRange, Graph, Neighbors)          
    npoints = np.divide(b, b + 1)*npoints + np.divide(count_pts, b + 1)
    Exp_while_loops = np.divide(b, b + 1)*Exp_while_loops + np.divide(count_while, b + 1)
    
    

    print(' Exp no of points :', npoints)
    print(' Expected while loops :', Exp_while_loops)


#### Conclusions
print('\n \n Parameters: \n Kappa_0 = ', Kappa_0, ', Gamma = ', gamma, '\nIntensity of PPP = ', Kappa, '\n Interaction Range = ', IntRange,'\n # iterations = ', Itot)
print('No of cells = ', K*K)
print('Results:')
print('  Estimated # points:', npoints)
print('  Expected # while loops:', Exp_while_loops)
print('Program Ending time: ', datetime.datetime.now().time())



     
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:55:29 2018

@author: uqsbabum
"""


# Discrete PRS alg to generate samples from Penetrable Spheres Mixture Model on [0,1]^2 with parameter kappa_1, Kappa_2 and interaction range,
# and abosolutely continuous with respect to
# the Poisson point process with intensity Kappa = (Kappa_1+Kappa_2)/(pi*IntRange*IntRange), where IntRange is
# the interaction range of the model. Note that type i points in PPP have intensity Kappa_i/(pi*IntRange*IntRange).

def reset():
    try:
        get_ipython().magic('reset -sf')  #analysis:ignore
    except NameError:
        pass
reset()


import numpy as np
import time
import datetime





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

###### Subroutine returns True if NC is within the Interaction range from an other type point in A. Used in dCFTP function.
def isBlocked(A, NC, IntRange):
    B = [x for x in A if x[2] != NC[2]]
    for x in B:
        if euclidean2dim(x[0:2], NC[0:2]) < IntRange:
            return True
    return False


###### Dominated CFTP method for generating a perfect sample of Strauss process on Cell numbered 'Cell_no' with the same parameters
def dCFTP(Cell_no, Kappa1, Kappa2, IntRange):
    x_cord = Cell_no%K
    y_cord = int(Cell_no/K)
    Cell_x_len = min(IntRange, 1 - x_cord*IntRange)
    Cell_y_len = min(IntRange, 1 - y_cord*IntRange)

    Kappa_cell = Cell_x_len*Cell_y_len*np.divide(4*(Kappa1 + Kappa2), np.pi*IntRange*IntRange)

    D = []
    events = [] #Records the events happened
    gamma = Kappa1/(Kappa1 + Kappa2)
    D.append([(Cell_x_len*np.random.random_sample(), Cell_y_len*np.random.random_sample(), np.random.binomial(1, gamma)) for i in range(np.random.poisson(Kappa_cell))]) # Intial state

    D.append(D[0][:])
    no_cirD = len(D[0]) #Number of disks at time 0

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
#                    print("Departure update of U and L")
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
#                    print("Arrival update of U and L")

                Block_in_U = isBlocked(U, events[n - j - 1][1], IntRange)
                Block_in_L = isBlocked(L, events[n - j - 1][1], IntRange)
#                    print((OL_ind_U, OL_ind_L))

                if Block_in_U == False:
                    L.append(events[n - j - 1][1])
                    U.append(events[n - j - 1][1])

                elif Block_in_L == False:
                    U.append(events[n - j - 1][1])

            if (set(L) <= set(U)) == False:
                print("L is not a subset of U")
                quit()
        ## Terminating the loop when U == L
        if len(L) == len(U):
            break


        ## Loop for doubling the length of the dominated process
        for j in range(n):
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
        State.append((pt[0] + x_cord*IntRange, pt[1] + y_cord*IntRange, pt[2]))

    return State


###### Node of the dependency graph
class GNode:
    def __init__(self, initdata):
        self.event = initdata
        self.next = None

###### Node for keeping a neighbor of a given node of the dependency graph
class NNode:
    def __init__(self, initdata):
        self.neighbor = initdata
        self.next = None


##### Subroutine for generating the nodes of the dependency graph
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


###### Subroutine for identifying the neighbors of each event in the dependency graph
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
def PSM_NPRS(K, Kappa1, Kappa2, IntRange, Graph, Neighbors):


    Disks1 = [[] for _ in range(K*K)]
    Disks2 = [[] for _ in range(K*K)]
    BadCells = [False for _ in range(K*K)]

    ## Intial state generation
    for x in range(K*K):
        State = dCFTP(x, Kappa1, Kappa2, IntRange)
        Disks1[x] = [x[0:2] for x in State if x[2] == 1]
        Disks2[x] = [x[0:2] for x in State if x[2] == 0]

    # Start of the while loop of Partial Rejection Sampling Method
    count_while = 0
    while True:
        count_while = count_while + 1

        ## Identifying the bad events
        stop = True
        pt = Graph

        while pt != None:
            (x,y) = list(pt.event)

            if Overlap(Disks1[x], Disks2[y], IntRange) or Overlap(Disks2[x], Disks1[y], IntRange):
                stop = False
                BadCells[x] = True
                BadCells[y] = True

            pt = pt.next

        ## PRS is terminated if there are no bad events
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
                        if (not ResCells[y]) and (len(Disks1[y]) > 0 or len(Disks2[y]) > 0):
                            stop = False
                            ResCells[y] = True
                        current = current.next



        for x in range(K*K):
            if (len(Disks1[x]) > 0 or len(Disks1[x]) > 0) and ResCells[x]:
                current = Neighbors[x]
                while current != None:
                    ResCells[current.neighbor] = True
                    current = current.next


        for x in range(K*K):
            if ResCells[x]:
                State = dCFTP(x, Kappa1, Kappa2, IntRange)
                Disks1[x] = [x[0:2] for x in State if x[2] == 1]
                Disks2[x] = [x[0:2] for x in State if x[2] == 0]

            BadCells[x] = False


    ## Conclusions after the iteration
    npoints1 = 0
    npoints2 = 0
    for x in range(K*K):
        npoints1 = npoints1 + len(Disks1[x])
        npoints2 = npoints2 + len(Disks2[x])

    return (npoints1, npoints2, count_while)



###### Model parameters
Kappa1 = 0.05
Kappa2 = 0.05
IntRange = 1/100

Itot = 1000  ### Number of samples generated

###### Partition parameters
K = int(1/IntRange)
Cell_edge = np.divide(1, K)
Kappa_0 = Kappa1 + Kappa2
Kappa = np.divide(4*Kappa_0, np.pi*IntRange*IntRange)
Kappa_cell = Kappa*Cell_edge*Cell_edge


###### Estimation parameters
Est_time = 0
Exp_while_loops = 0
Exp_total_pts = 0
npoints1 = 0
npoints2 = 0



print('\n \n Parameters: \n Kappa1 = ', Kappa1, ', Kappa2 = ', Kappa2, '\n Intensity of PPP type 1 =', np.divide(4*Kappa1, np.pi*IntRange*IntRange), '\n Intensity of PPP type 2 = ', np.divide(4*Kappa2, np.pi*IntRange*IntRange), '\n Total Intensity = ', Kappa, '\n Interaction Range = ', IntRange,'\n # iterations = ', Itot)
print('No of cells = ', K*K, '\n Intensity per cell of PPP: ', Kappa_cell)

start_time = time.time()
np.random.seed(0)
print('Program Starting time: ', datetime.datetime.now().time())
print('++++++++++ New PRS Alg for Strauss Process +++++++++++++++')


Graph = generateVertices(K)
Neighbors = generateNeighbors(K)

for b in range(Itot):
#    print('---------- New PRS Iteration: ', b+1, '---------')


    (count_pts1, count_pts2, count_while) = PSM_NPRS(K, Kappa1, Kappa2, IntRange, Graph, Neighbors)
    npoints1 = np.divide(b, b + 1)*npoints1 + np.divide(count_pts1, b + 1)
    npoints2 = np.divide(b, b + 1)*npoints2 + np.divide(count_pts2, b + 1)
    Exp_while_loops = np.divide(b, b + 1)*Exp_while_loops + np.divide(count_while, b + 1)


#
#    print(' Exp no of type 1 points :', npoints1)
#    print(' Exp no of type 2 points :', npoints2)
#    print(' Expected while loops :', Exp_while_loops)

end_time = time.time()
Est_time = end_time - start_time

print('\n \n Parameters: \n Kappa1 = ', Kappa1, ', Kappa2 = ', Kappa2, '\n Intensity of PPP type 1 =', np.divide(4*Kappa1, np.pi*IntRange*IntRange), '\n Intensity of PPP type 2 = ', np.divide(4*Kappa2, np.pi*IntRange*IntRange), '\n Total Intensity = ', Kappa, '\n Interaction Range = ', IntRange,'\n # iterations = ', Itot)
print('No of cells = ', K*K, '\n Intensity per cell of PPP: ', Kappa_cell)
print('Results:')
print(' Exp no of type 1 points :', npoints1)
print(' Exp no of type 2 points :', npoints2)
print('  Expected # while loops:', Exp_while_loops)
#print('Expected total # points :', Exp_total_pts)
print('Time Taken: ', Est_time)
##plt.plot(Beta, Est_time)
##plt.xlabel('Intensity of the Poisson point process')
##plt.ylabel('Time taken to generate 10,000 samples (in secs)')
##plt.show()
#

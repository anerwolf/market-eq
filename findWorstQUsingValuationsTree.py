from treeNode import treeNode
from utils import *
from scipy.optimize import linprog
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

Path('./logs-dir').mkdir(parents=True, exist_ok=True)
logsDir = './logs-dir/'
nonFeasibleMsgCnt = 500

mMax = 100
m = 3
eps = 0.000001
binarySearchPercent = 1.01

worstQ = [0] * m

cStart = m + eps

while m <= mMax:
    Path(logsDir + 'm_' + str(m) + '_figs').mkdir(parents=True, exist_ok=True)
    f = open(logsDir + 'm_' + str(m) + '_log.txt', "w")
    
    stoppingThreshold = 0.0005
    foundWorstQ = False
    print('Searching for worst Q. m: ' + str(m))
    cMin = 0
    cMax = m + eps
    
    plotBool = True

    #Build tree
    root = treeNode(0, m, 1, 1, None, [], [], [], [])
    root.createChildren()
    currentChildren = root.children
    for k in range(1, m + 1):
        numOfChildren = len(currentChildren)
        nextChildren = []
        for childIdx in range(numOfChildren):
            currentChildren[childIdx].createChildren()
            nextChildren = nextChildren + currentChildren[childIdx].children
        currentChildren = nextChildren
    leafs = currentChildren
    minLeafIdx = 0
    print('Number of leafs: ' + str(len(leafs)))
    f.write('Number of leafs: ' + str(len(leafs))+ '\n')

    cnt = -1
    while not foundWorstQ:
        cnt += 1

        if cnt == 0:
            c = cStart
        else:
            cTmp = ((cMax - cMin) / 2) + cMin
            c = min(cTmp, binarySearchPercent * c)

        print('Current c: ' + str(c) + '. cMin: ' + str(cMin) + '. cMax: ' + str(cMax))
        f.write('Current c: ' + str(c) + '. cMin: ' + str(cMin) + '. cMax: ' + str(cMax)+ '\n')
        

        #Solve Lp`s
        obj = list([0] * (2 * (m + 1)))
        bnd = list([((0, float("inf")))]  * (2 * (m + 1)))
        foundFeasibleSol = False
        
        for leafIdx in range(minLeafIdx, len(leafs)):
            LHSConstraints, RHSConstraints, LHSConstraintsSym, RHSConstraintsSym, LHSEqConstraints, RHSEqConstraints = leafs[leafIdx].getConstraints()
            LHSConstraintsWithC, RHSConstraintsWithC = addCToConstraints(LHSConstraintsSym, RHSConstraintsSym, c, m)
            LHSConstraints = LHSConstraints + LHSConstraintsWithC
            RHSConstraints = RHSConstraints + RHSConstraintsWithC
            try:
                opt = linprog(c=obj, A_ub=LHSConstraints, b_ub=RHSConstraints, A_eq=LHSEqConstraints, b_eq=RHSEqConstraints, bounds=bnd)
            except:
                if leafIdx == minLeafIdx or leafIdx % nonFeasibleMsgCnt == 0:
                    print('Idx: ' + str(leafIdx) + '. No feasible solution for current leaf. c: ' + str(c) + ', m: ' + str(m))
            else:
                if opt.status == 0 and opt.success:
                    foundFeasibleSol = True
                    v1 = opt.x[0:m + 1]
                    v2 = opt.x[m+1:]
                    f1, f2, f1PlusF2, sw, q = calculateF1F2SWQ(v1, v2, m)
                    maxSWIdx = np.argmax(sw)
                    print('Idx: ' + str(leafIdx) + '. Feasible solution found.  c: ' + str(c) + ', m: ' + str(m))
                    f.write('Idx: ' + str(leafIdx) + '. Feasible solution found.  c: ' + str(c) + ', m: ' + str(m)+ '\n')
                    print('v1: ' + str(v1))
                    f.write('v1: ' + str(v1)+ '\n')
                    print('v2: ' + str(v2))
                    f.write('v2: ' + str(v2)+ '\n')
                    print('sw: ' + str(sw))
                    f.write('sw: ' + str(sw)+ '\n')
                    print('q: ' + str(q))
                    f.write('q: ' + str(q)+ '\n')
                    print('sw max idx (k) : ' + str(maxSWIdx) + ', max sw: ' + str(sw[maxSWIdx]) + ', q: ' + str(q[maxSWIdx]))
                    f.write('sw max idx (k) : ' + str(maxSWIdx) + ', max sw: ' + str(sw[maxSWIdx]) + ', q: ' + str(q[maxSWIdx])+ '\n')
                    print('f1: ' + str(f1))
                    f.write('f1: ' + str(f1)+ '\n')
                    print('f2: ' + str(f2))
                    f.write('f2: ' + str(f2)+ '\n')
                    print('f1 + f2: ' + str(f1PlusF2))
                    f.write('f1 + f2: ' + str(f1PlusF2)+ '\n')
                    if plotBool:
                        plt.plot(sw, label='v1[m-k] + v2[k]')
                        plt.plot(f1PlusF2, label='f1[k] + f2[k]')
                        plt.xlabel('k')
                        plt.legend()
                        plt.savefig(logsDir + 'm_' + str(m) + '_figs/' + 'm_' + str(m) + '_c_' + str(c) + 'VAndF.png')
                        plt.close()
                        v1Add = turnToAdditive(v1)
                        v2Add = turnToAdditive(v2)
                        plt.plot(v1, label='v1')
                        plt.plot(v1Add, label='v1 Additive bound')
                        plt.plot(v2, label='v2')
                        plt.plot(v2Add, label='v2 Additive bound')
                        plt.xlabel('k')
                        plt.legend()
                        plt.savefig(logsDir + 'm_' + str(m) + '_figs/' + 'm_' + str(m) + '_c_' + str(c) + 'V1V2.png')
                        plt.close()
                    break

        if foundFeasibleSol:
            minLeafIdx = leafIdx
            if (c - cMin < stoppingThreshold) and (cMax - c < stoppingThreshold):
                print('Found worst Q: ' + str(c))
                f.write('Found worst Q: ' + str(c)+ '\n')
                f.close()
                worstQ.append(c)
                cStart = c
                m += 1
                foundWorstQ = True
            else: 
                cMin = c
        else:
            cMax = c
        
    ones = [1] * len(worstQ)
    plt.plot(worstQ, label='q')
    plt.plot(ones, color='r')
    plt.xlabel('m')
    plt.legend()
    plt.savefig(logsDir + 'maxM_' + str(m-1) + '_worstQFig.png')
    plt.close()
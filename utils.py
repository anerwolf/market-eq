def calculateF1F2SWQ(v1, v2, m):
    f1 = [0] * len(v1)
    f2 = [0] * len(v2)
    f1PlusF2 = [0] * len(v2)
    sw = [0] * len(v2)
    q = [0] * len(v2)
    for k in range(m + 1):
        sw[k] = v1[m - k] + v2[k]
        for l in range(k):
            f_1_tmp = k * (v1[l + 1 + m - k] - v1[m - k]) / (l + 1)
            if f_1_tmp > f1[k]:
                f1[k] = f_1_tmp
        for l in range(m - k):
            f_2_tmp = (m - k) * (v2[l + 1 + k] - v2[k]) / (l + 1)
            if f_2_tmp > f2[k]:
                f2[k] = f_2_tmp
        f1PlusF2[k] = f1[k] + f2[k]
        q[k] = f1PlusF2[k] / sw[k]
    return f1, f2, f1PlusF2, sw, q

def addCToConstraints(LHSConstraintsSym, RHSConstraintsSym, c, m):
    LHSConstraintsWithC = []
    RHSConstraintsWithC = []
    for constraintIdx in range(len(LHSConstraintsSym)):
        constraint = list([0] * (2 * (m + 1)))
        for i in range(len(LHSConstraintsSym[constraintIdx])):
            idx = LHSConstraintsSym[constraintIdx][i][0]
            string = LHSConstraintsSym[constraintIdx][i][1]
            value = LHSConstraintsSym[constraintIdx][i][2]
            constraint[idx] = value
            if string == 'c':
                constraint[idx] += c
        LHSConstraintsWithC.append(constraint)
        RHSConstraintsWithC.append(RHSConstraintsSym[constraintIdx])
    
    return LHSConstraintsWithC, RHSConstraintsWithC
class treeNode:
    
    def __init__(self, k, m, l1Cap, l2Cap, parent, LHSConstraints, RHSConstraints, LHSConstraintsSym, RHSConstraintsSym):
        self.k = k
        self.m = m
        self.l1Cap = l1Cap
        self.l2Cap = l2Cap
        self.parent = parent
        self.children = []
        self.isLeaf = k == m + 1
        self.isRoot = k == 0
        self.LHSConstraints = LHSConstraints
        self.RHSConstraints = RHSConstraints
        self.LHSConstraintsSym = LHSConstraintsSym
        self.RHSConstraintsSym = RHSConstraintsSym
        self.eps = 0.000001
        if (self.isLeaf):
            self.createSAConstarints()
            self.createMonotoneConstraints()

    def createChildren(self):
        start1 = 1
        if (self.k == 0):
            l1Zero = True
            end1 = 1
        else: 
            l1Zero = False
            end1 = self.k
        start2 = 1
        if (self.m - self.k == 0):
            l2Zero = True
            end2 = 1
        else: 
            l2Zero = False
            end2 = self.m - self.k

        if (self.m - self.k < self.l1Cap):
            start1 = self.l1Cap - (self.m - self.k)
            end1 = self.l1Cap - (self.m - self.k)
        if (self.k < self.l2Cap):
            start2 = self.l2Cap - self.k
            end2 = self.l2Cap - self.k
            
        v2StartIdx = self.m + 1
        
        for l1 in range(start1, end1 + 1):
            for l2 in range(start2, end2 + 1):
                LHSConstraintsSym = list(self.LHSConstraintsSym)
                RHSConstraintsSym = list(self.RHSConstraintsSym)
                constraintSym = []
                if l1Zero:
                    # constraint[self.m - self.k] = self.c
                    constraintSym.append([self.m - self.k, 'c', 0])
                else: 
                    # constraint[self.m - self.k] = self.c + (self.k / l1)
                    constraintSym.append([self.m - self.k, 'c', (self.k / l1)])
                    # constraint[l1 + self.m - self.k] = -self.k / l1
                    constraintSym.append([l1 + self.m - self.k, '', -self.k / l1])
                if l2Zero:
                    # constraint[v2StartIdx + self.k] = self.c
                    constraintSym.append([v2StartIdx + self.k, 'c', 0])
                else:
                    # constraint[v2StartIdx + self.k] = self.c + ((self.m - self.k) / l2)
                    constraintSym.append([v2StartIdx + self.k, 'c', ((self.m - self.k) / l2)])
                    # constraint[v2StartIdx + l2 + self.k] = -(self.m - self.k) / l2
                    constraintSym.append([v2StartIdx + l2 + self.k, '', -(self.m - self.k) / l2])
                LHSConstraintsSym.append(constraintSym)
                RHSConstraintsSym.append(-self.eps)

                LHSConstraints = list(self.LHSConstraints)
                RHSConstraints = list(self.RHSConstraints)
                for j1 in range(1, self.k + 1):
                    if (l1 == j1):
                        continue
                    constraint = list([0] * (2 * (self.m + 1)))
                    constraint[j1 + self.m - self.k] = 1 / j1
                    constraint[self.m - self.k] = (1 / l1) - (1 / j1)
                    constraint[l1 + self.m - self.k] = -1 / l1
                    LHSConstraints.append(constraint)
                    RHSConstraints.append(0)
                for j2 in range(1, self.m - self.k + 1):
                    if (l2 == j2):
                        continue
                    constraint = list([0] * (2 * (self.m + 1)))
                    constraint[v2StartIdx + j2 + self.k] = 1 / j2
                    constraint[v2StartIdx + self.k] = (1 / l2) - (1 / j2)
                    constraint[v2StartIdx + l2 + self.k] = -1 / l2
                    LHSConstraints.append(constraint)
                    RHSConstraints.append(0)
                self.children.append(treeNode(self.k + 1, self.m, l1, l2, self, LHSConstraints, RHSConstraints, LHSConstraintsSym, RHSConstraintsSym))

    def createSAConstarints(self):
        v2StartIdx = self.m + 1
        for i in range(1, self.m):
            for j in range(1, min(i + 1, self.m - i + 1)):
                constraint = list([0] * (2 * (self.m + 1)))
                constraint[i + j] = 1
                if i == j:
                    constraint[i] = -2
                else:    
                    constraint[i] = -1
                    constraint[j] = -1
                self.LHSConstraints.append(constraint)
                self.RHSConstraints.append(0)
                constraint = list([0] * (2 * (self.m + 1)))
                constraint[v2StartIdx + i + j] = 1
                if i == j:
                    constraint[v2StartIdx + i] = -2
                else:
                    constraint[v2StartIdx + i] = -1
                    constraint[v2StartIdx + j] = -1
                self.LHSConstraints.append(constraint)
                self.RHSConstraints.append(0)

    def createMonotoneConstraints(self):
        v2StartIdx = self.m + 1
        for i in range(1, self.m):
            constraint = list([0] * (2 * (self.m + 1)))
            constraint[i] = 1
            constraint[i + 1] = -1
            self.LHSConstraints.append(constraint)
            self.RHSConstraints.append(0)
            constraint = list([0] * (2 * (self.m + 1)))
            constraint[v2StartIdx + i] = 1
            constraint[v2StartIdx + i + 1] = -1
            self.LHSConstraints.append(constraint)
            self.RHSConstraints.append(0)

    def getConstraints(self):
        v2StartIdx = self.m + 1
        #Create norm constraints
        LHSEqConstraints = []
        RHSEqConstraints = []
        constraint = list([0] * (2 * (self.m + 1)))
        constraint[0] = 1
        LHSEqConstraints.append(constraint)
        RHSEqConstraints.append(0)
        constraint = list([0] * (2 * (self.m + 1)))
        constraint[v2StartIdx] = 1
        LHSEqConstraints.append(constraint)
        RHSEqConstraints.append(0)
        return self.LHSConstraints, self.RHSConstraints, self.LHSConstraintsSym, self.RHSConstraintsSym, LHSEqConstraints, RHSEqConstraints

                    
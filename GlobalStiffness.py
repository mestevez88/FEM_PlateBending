import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

class FEM_Assemble():
    def __init__(self, VoidCheck, LX, LY, EX, EY, LC_Nodes, BC, Gpx, Gpy, t):
        self.VoidCheck = VoidCheck
        self.LX = LX
        self.LY = LY
        self.EX = EX
        self.EY = EY
        self.DGoF = 3*(EX+1)*(EY+1)
        self.Elements = EX*EY
        self.LC_Nodes = LC_Nodes
        self.BC = BC
        self.Gpx = Gpx
        self.Gpy = Gpy
        self.t = t
    def Global_Stiffnes_Matrix(self, ElementNodes, Kij):
        'Number of Degrees of Freedom of global stiffness matrix'
        rows, cols, data = [], [], []
        'Loop through all elementes e'
        for e in range(self.Elements):
            'Get the node index for the element e'
            idx = ElementNodes[e].astype(int)
            'Element DGoF [u1,v1,w1, u2,v2,w2, u3,v3,w3, u4,v4,w4]'
            elementDof = np.empty(12, dtype=int)
            'Get the global degrees of freedom for the element nodes'
            elementDof[0::3] = 3*idx
            elementDof[1::3] = 3*idx + 1
            elementDof[2::3] = 3*idx + 2
            'Create the local stiffness matrix k_loc (12x12) for element e'
            k_loc = np.zeros((12, 12), dtype=np.float64)
            'Calculate the value of integral (3) by adding the terms running through indicies (i,j)'
            for i in range(self.Gpx):
                for j in range(self.Gpy):
                    if self.VoidCheck[e] == 1:
                        k_loc += Kij[e, i, j, 1, :, :] #solid
                    else:
                        k_loc += Kij[e, i, j, 0, :, :] #voided
            'Create meshgrid with Degrees of Freedom'
            rr, cc = np.meshgrid(elementDof, elementDof, indexing='ij')
            'Flatten the meshgrid indicies'
            rows.append(rr.ravel())
            cols.append(cc.ravel())
            'Flatten k_loc entries'
            data.append(k_loc.ravel())
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        data = np.concatenate(data)
        'Global stiffnes matrix K (DGoF, DGoF)'
        'Assembly of the global stiffness matrix K by inserting the data in the specified matrix entry (rows, cols)'
        K = coo_matrix((data, (rows, cols)), shape=(self.DGoF, self.DGoF)).tocsc()
        return K
    def Displacement_Vector(self, K, F):
        # Build active dof mask directly (more efficient than list operations)
        activeDof_mask = np.ones(self.DGoF, dtype=bool)
        activeDof_mask[self.BC] = False   
        # Extract active block using boolean indexing
        ActiveStiffness = K[activeDof_mask][:, activeDof_mask]
        # Build active force directly using mask
        ActiveForce = F[activeDof_mask, 0]  # Extract as 1D array directly
        displacements = np.zeros((self.DGoF, 1), dtype=np.float64)
        U_active = spsolve(ActiveStiffness, ActiveForce)  # (n_active,)
        displacements[activeDof_mask, 0] = U_active
        return displacements

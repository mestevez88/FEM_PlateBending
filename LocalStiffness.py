import numpy as np

class FEM_PlateBending():
    def __init__(self, LX, LY, EX, EY, Gpx, Gpy, t, materials):
        self.LX = LX
        self.LY = LY
        self.EX = EX
        self.EY = EY
        self.Elements = EX*EY
        self.Gpx = Gpx
        self.Gpy = Gpy
        self.Elements = EX*EY
        self.t = t
        self.materials = materials
        'Location of Gauss Points: The order in the double for loop for Gauss points (xi, eta) is the following: (-,-), (-,+), (+,-), (+,+)'
        'This means that the noode coordinates of the mesh are in the following order: LeftBottom, LeftUp, RightBottom, RightUp'
        'Location of Gauss Points in the X (Xi) direction'
        self.Xi={1:[0], 2:[-np.sqrt(1/3), np.sqrt(1/3)]}
        'Gauss Weights on the X direction'
        self.WXi={1:[2], 2:[1, 1]}
        'Location of Gauss Points in the Y (Eta) direction'
        self.Eta={1:[0],2:[-np.sqrt(1/3), np.sqrt(1/3)]}
        'Gauss Weights on the Y direction'
        self.WEta={1:[2], 2:[1, 1]}
        'Materials'
        E_sol, nu_sol = float(self.materials[0,0]), float(self.materials[0,1])
        E_void, nu_void = float(self.materials[1,0]), float(self.materials[1,1])
        def Dmat(E, nu):
            return (E*(self.t**3) / (12*(1.0 - nu**2))) * np.array(
                [[1.0, nu, 0.0],
                [nu, 1.0, 0.0],
                [0.0, 0.0, (1.0 - nu) / 2.0]],
                dtype=np.float64
            )
        self.D_sol  = Dmat(E_sol,  nu_sol)
        self.D_void = Dmat(E_void, nu_void)
    'Derivative of the Shape functions'
    def Nder_ij(self, xi, eta):
        N_ij = (1/4)*np.array([[-(1-eta), -(1+eta), 1-eta, 1+eta], [-(1-xi), 1-xi, -(1+xi), 1+xi]], dtype=np.float64)
        return(N_ij)
    'Shape function Matrix'
    def Nij(self, xi, eta):
        Nij = (1/4)*np.array([(1-xi)*(1-eta), (1-xi)*(1+eta), (1+xi)*(1-eta), (1+xi)*(1+eta)], dtype=np.float64)
        return(Nij)
    'Jacobian determinant'
    def Jacobian(self, nodeCoor, Nder):
        J_ij = np.matmul(Nder, nodeCoor)
        Jinv_ij = np.linalg.inv(J_ij)
        XYDerivative = np.matmul(Jinv_ij, Nder)
        detJ_ij = np.linalg.det(J_ij)
        return(J_ij, Jinv_ij, XYDerivative, detJ_ij)
    'Element Stiffness Matrices'
    def Local_Stiffness(self, ElementNodes, NodeCoord):
        lenGpx = len(self.Xi[self.Gpx])
        lenGpy = len(self.Eta[self.Gpy])
        WiWj = np.zeros((lenGpx, lenGpy), dtype=np.float64)
        for i in range(lenGpx):
            Wi = self.WXi[self.Gpx][i]
            for j in range(lenGpy):
                Wj = self.WEta[self.Gpy][j]
                WiWj[i, j] = Wi * Wj
        # Allocate outputs
        Kij = np.zeros((self.Elements, lenGpx, lenGpy, 2, 12, 12), dtype=np.float64)
        Bij = np.zeros((self.Elements, lenGpx, lenGpy, 3, 12), dtype=np.float64)
        # Loop over elements
        for e in range(self.Elements):
            idx = ElementNodes[e].astype(int)
            Q = NodeCoord[idx, :]
            #centroidal coordinates
            x_cent = np.mean(Q[:,0])
            y_cent = np.mean(Q[:,1])
            el_xy_cent = np.array([[x_cent, y_cent],[x_cent, y_cent],[x_cent, y_cent],[x_cent, y_cent]], dtype=np.float64)
            #Build C matrix
            C = np.zeros((12,12), dtype=np.float64)
            'Perform double for loop in correspondence with Gauss ponts'
            for i in range(lenGpx):
                xi = self.Xi[self.Gpx][i]
                for j in range(lenGpy):
                    eta = self.Eta[self.Gpy][j]
                    #centroidal coordinates
                    x = Q[2*i+j,0] - x_cent
                    y = Q[2*i+j,1] - y_cent
                    index = 3*(2*i + j)
                    C[index, 0] = 1 ; C[index, 1] = x ; C[index, 2] = y ; C[index, 3] = x**2 ; C[index, 4] = x*y ; C[index, 5] = y**2 ; C[index, 6] = x**3 ; C[index, 7] = x**2*y ; C[index, 8] = x*y**2 ; C[index, 9] = y**3 ; C[index, 10] = x**3*y ; C[index, 11] = x*y**3
                    C[index + 1, 2] = 1 ; C[index + 1, 4] = x ; C[index + 1, 5] = 2*y ; C[index + 1, 7] = x**2 ; C[index + 1, 8] = 2*x*y ; C[index + 1, 9] = 3*y**2 ; C[index + 1, 10] = x**3 ; C[index + 1, 11] = 3*x*y**2
                    C[index + 2, 1] = -1 ; C[index + 2, 3] = -2*x; C[index + 2, 4] = -y ; C[index + 2, 6] = -3*x**2; C[index + 2, 7] = -2*x*y; C[index + 2, 8] = -y**2 ; C[index + 2, 10] = -3*x**2*y; C[index + 2, 11] = -y**3
            # Precompute B for each Gauss point using cached Jinv
            for i in range(lenGpx):
                xi = self.Xi[self.Gpx][i]
                for j in range(lenGpy):
                    eta = self.Eta[self.Gpy][j]
                    'N_ij derivatives'
                    Nder = self.Nder_ij(xi,eta)
                    J = Nder @ (Q - el_xy_cent) #Computing J after translating the centroid of the plate to the origin.
                    detJ = float(np.linalg.det(J))
                    #centroidal coordinates
                    x_int = Q[:,0] - x_cent
                    y_int = Q[:,1] - y_cent
                    'Real coordinates of Gauss Points, this simply multiplies the GaussLoc in Xi (resp. Eta) by the corresponding mesh coordinate X (resp. Y)'
                    x = self.Nij(xi, eta) @ x_int
                    y = self.Nij(xi, eta) @ y_int
                    #Build H for each node (3x12)
                    H = np.zeros((3, 12), dtype=np.float64)
                    H[0,3] = 2; H[0,6] = 6*x; H[0,7] = 2*y; H[0,10] = 6*x*y
                    H[1,5] = 2; H[1,8] = 2*x; H[1,9] = 6*y; H[1,11] = 6*x*y
                    H[2,4] = 2; H[2,7] = 4*x; H[2,8] = 4*y; H[2,10] = 6*x**2; H[2,11] = 6*y**2
                    C_inv = np.linalg.inv(C)
                    # Build B 
                    B = H @ C_inv
                    Bij[e, i, j, :, :] = B
                    w = WiWj[i, j] * detJ
                    BT = B.T
                    # Void & solid element matrices for this GP
                    Kij[e, i, j, 0, :, :] = w * (BT @ (self.D_void @ B))
                    Kij[e, i, j, 1, :, :] = w * (BT @ (self.D_sol  @ B))
        return Kij

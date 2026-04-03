import numpy as np
from plots import Plot_PlateCenter

class Post_Processing():
    def __init__(self, VoidCheck, LX, LY, EX, EY, Gpx, Gpy, t, materials, Bij, NodeCoord, ElementNodes, v):
        self.VoidCheck = VoidCheck
        self.LX = LX
        self.LY = LY
        self.EX = EX
        self.EY = EY
        self.Elements = EX*EY
        self.dx = self.LX/self.EX
        self.dy = self.LY/self.EY
        self.Gpx = Gpx
        self.Gpy = Gpy
        self.t = t
        self.materials = materials
        self.Bij = Bij
        self.v = v
        self.NodeCoord = NodeCoord
        self.ElementNodes = ElementNodes
        'Location of Gauss Points in the X (Xi) direction'
        self.Xi={1:[0], 2:[-np.sqrt(1/3), np.sqrt(1/3)]}
        'Gauss Weights on the X direction'
        self.WXi={1:[2], 2:[1, 1]}
        'Location of Gauss Points in the Y (Eta) direction'
        self.Eta={1:[0],2:[-np.sqrt(1/3), np.sqrt(1/3)]}
        'Gauss Weights on the Y direction'
        self.WEta={1:[2], 2:[1, 1]}
        "Materials"
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
        self.x = NodeCoord[:,0]
        self.y = NodeCoord[:,1]
    def Calculate_Forces(self, plot):
        Curv = np.zeros((self.Elements, self.Gpx, self.Gpy, 3, 1), dtype=np.float64)
        M = np.zeros((self.Elements, self.Gpx, self.Gpy, 3, 1), dtype=np.float64)
        MX = np.zeros(self.Elements, dtype=np.float64)
        MY = np.zeros(self.Elements, dtype=np.float64)
        MXY = np.zeros(self.Elements, dtype=np.float64)
        for e in range(self.Elements):
            idx = self.ElementNodes[e].astype(int)
            idx = 3*idx
            idx_ = np.concatenate((idx, idx+1))
            idx_ = np.concatenate((idx_, idx+2))
            idx_.sort()
            for i in range(self.Gpx):
                for j in range(self.Gpy):
                    kappa = self.Bij[e, i, j, :, :] @ self.v[idx_]
                    Curv[e, i, j, :, :] = kappa
                    if self.VoidCheck[e] == 1:
                        M[e, i, j, :, :] = self.D_sol @ kappa
                    else:
                        M[e, i, j, :, :] = self.D_void @ kappa
            'Take the mean over the all Gauss Points'
            MX[e] = .25*(M[e, 0, 0, 0, 0] + M[e, 0, 1, 0, 0] + M[e, 1, 0, 0, 0] + M[e, 1, 1, 0, 0])
            MY[e] = .25*(M[e, 0, 0, 1, 0] + M[e, 0, 1, 1, 0] + M[e, 1, 0, 1, 0] + M[e, 1, 1, 1, 0])
            MXY[e] = .25*(M[e, 0, 0, 2, 0] + M[e, 0, 1, 2, 0] + M[e, 1, 0, 2, 0] + M[e, 1, 1, 2, 0])
        MXY_2d = np.reshape(MXY, (self.EX, self.EY))
        'Shear in X'
        VX = np.zeros(self.Elements, dtype=np.float64)
        MX_2d = np.reshape(MX, (self.EX, self.EY))
        VX = np.gradient(MX_2d, self.dx, axis=0) + np.gradient(MXY_2d, self.dy, axis=1)
        'Shear in Y'
        VY = np.zeros(self.Elements, dtype=np.float64)
        MY_2d = np.reshape(MY, (self.EX, self.EY))
        VY = np.gradient(MY_2d, self.dy, axis=1) + np.gradient(MXY_2d, self.dx, axis=0)
        if plot:
            Plot_PlateCenter(self.x, self.y, self.EX, self.EY, self.LX, self.LY, MX, 'Moment in X-direction (MX, kgf*cm/cm)')
            Plot_PlateCenter(self.x, self.y, self.EX, self.EY, self.LX, self.LY, MY, 'Moment in Y-direction (MY, kgf*cm/cm)')
            Plot_PlateCenter(self.x, self.y, self.EX, self.EY, self.LX, self.LY, MXY, 'Torsion (MXY, kgf*cm/cm)')
            Plot_PlateCenter(self.x, self.y, self.EX, self.EY, self.LX, self.LY, VX, 'Shear X (VX, kgf)')
            Plot_PlateCenter(self.x, self.y, self.EX, self.EY, self.LX, self.LY, VY, 'Shear Y (VY, kgf)')





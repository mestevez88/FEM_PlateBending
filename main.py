import numpy as np
from LocalStiffness import FEM_PlateBending
from GlobalStiffness import FEM_Assemble
from PostProcessing import Post_Processing
from meshing import rectangularmesh
from Utilities import CCCC, SSSS, V_Beam, H_Beam, UniformLoad
from plots import Plot_PlateNodes

'Geometric Properties'
LX = 1500
LY = 1500
EX = 3*15
EY = 3*15
Gpx = 2
Gpy = 2
t = 3

'Create a rectangular mesh of sides LX, LY with EX and EY subdivisions with Lower Left mesh node at coordinate (0,0)'
RectOutput = rectangularmesh(LX, LY, EX, EY, 0, 0)
NodeCoord = RectOutput[0]
ElementNodes = RectOutput[1]
'Extract Node x and y coordinates'
x = NodeCoord[:,0]
y = NodeCoord[:,1]

'Introduce Material Properties'
materials = np.zeros((2,2), dtype=np.float64)
'Solid Material Properties'
materials[0,0] = 2_039_000 #Modulus of Elasticity Steel (kgf/cm2)
materials[0,1] = 0.3 #Poisson's Ratio
'Void Material'
materials[1,0] = 1 #Modulus of Elasticity (1e-6 MPa)
materials[1,1] = 1e-6 #Poisson's Ratio

'Calculate the Element Stiffness Matrices Kij and Shape Matrices Bij'
K_Local = FEM_PlateBending(LX, LY, EX, EY, Gpx, Gpy, t, materials)
Kij, Bij = K_Local.Local_Stiffness(ElementNodes, NodeCoord)

'Specify Voided or Solid Elements'
VoidCheck = [1]*EX*EY #All solid elements

'Uniform Load (kg/cm2)'
LC_Nodes = UniformLoad(EX, EY, LX, LY, .01) #q (kg/cm**2), LX,LY (cm)
'Input Uniform Load values into Force Vector F'
DGoF = 3*(EX+1)*(EY+1) #Degrees of Freedom
F = np.zeros((DGoF, 1), dtype=np.float64)
for i in LC_Nodes:
    F[3*i] = LC_Nodes[i][0]

'Plate Boundary Conditions'
'LB(SS), LU(SS), RU(SS)'
# BC = [0, 1, 2, 3*EY, 3*EY+1, 3*EY + 2]
'Simply Supported on its four edges (SSSS)'
N_SSSS, BC_SSSS = SSSS(EX, EY) #Simply supported on its four edges.
'Clamped on its four edges (CCCC)'
# N_CCCC, BC_CCCC = CCCC(EX, EY) #Clamped on its four edges.
'Vertical Beam'
_, BC_V1 = V_Beam(V_Axes=[int(EX/3)], EY = EY, Rows = "All" , BC = "Simply Supported")
_, BC_V2 = V_Beam(V_Axes=[int(2*EX/3)], EY = EY, Rows = "All" , BC = "Simply Supported")
'Horizontal Beam'
_, BC_H1 = H_Beam(H_Axes=[int(EY/3)], EX = EX, Cols = "All", BC = "Simply Supported")
_, BC_H2 = H_Beam(H_Axes=[int(2*EY/3)], EX = EX, Cols = "All", BC = "Simply Supported")
'Combine Boundary Conditions'
# args = [BC_CCCC, BC_V1, BC_H1, BC_V2, BC_H2] #Clampled perimeter and intermediate beams
args = [BC_SSSS, BC_V1, BC_H1, BC_V2, BC_H2] #Simply supported perimeter and intermediate beams
# args = [BC_V1, BC_H1, BC_V2, BC_H2] #Free perimeter and intermediate beams
'Concatenate list'
BC  = np.concatenate(args, axis=0)
'Eliminate duplicate indicies'
BC = np.sort(np.unique(BC))

'Calculate vertical displacement vector'
K_Global = FEM_Assemble(VoidCheck, EX, EY, Gpx, Gpy, BC)
K = K_Global.Global_Stiffnes_Matrix(ElementNodes, Kij)
v = K_Global.Displacement_Vector(K,F)
'Plotting Displacements'
Plot_PlateNodes(x, y, EX, EY, LX, LY, v, 'Deflection (w, cm)')
PostProcessing = Post_Processing(VoidCheck, LX, LY, EX, EY, Gpx, Gpy, t, materials, Bij, NodeCoord, ElementNodes, v)
PostProcessing.Calculate_Forces(plot=True)
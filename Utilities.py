import numpy as np

'Boundary Conditions'
def CCCC(EX, EY): #Left, Right, Top, Bottom
    Nodes = []
    NX = (EX+1)
    NY = (EY+1)
    for i in range(NX*NY):
        if np.floor(i/NY) == 0 or np.floor(i/NY) == NX-1 or i%NY == 0 or i%NY == NY - 1:
            Nodes.append(i)
    Nodes = np.array(Nodes, dtype=int)
    DGoF = np.concatenate((3*Nodes, 3*Nodes+1, 3*Nodes+2), axis = None)
    DGoF = np.sort(DGoF)
    return(Nodes, DGoF)

def SSSS(EX, EY): #Left, Right, Top, Bottom
    Nodes = []
    NX = (EX+1)
    NY = (EY+1)
    for i in range(NX*NY):
        if np.floor(i/NY) == 0 or np.floor(i/NY) == NX-1 or i%NY == 0 or i%NY == NY-1:
            Nodes.append(i)
    Nodes = np.array(Nodes, dtype=int)
    DGoF = 3*Nodes
    DGoF = np.sort(DGoF)
    return(Nodes, DGoF)

def V_Beam(V_Axes, EY, Rows, BC):
    Nodes = []
    NY = (EY+1)
    if BC == "Simply Supported":
        for x in V_Axes:
            if Rows == "All":
                for y in range(NY):
                    v = x*NY + y
                    Nodes.append(v)
            else:
                for y in Rows:
                    v = x*NY + y
                    Nodes.append(v)
    Nodes = np.array(Nodes, dtype=int)
    DGoF = 3*Nodes
    DGoF = np.sort(DGoF)
    return(Nodes, DGoF)

def H_Beam(H_Axes, EX, Cols, BC):
    Nodes = []
    NX = (EX+1)
    if BC == "Simply Supported":
        for y in H_Axes:
            if Cols == "All":
                for x in range(NX):
                    v = y + x*NX
                    Nodes.append(v)
            else:
                for x in Cols:
                    v = y + x*NX
                    Nodes.append(v)
    Nodes = np.array(Nodes, dtype=int)
    DGoF = 3*Nodes
    DGoF = np.sort(DGoF)
    return(Nodes, DGoF)

'Load Functions'
def UniformLoad(EX, EY, LX, LY, q): #Left, Right, Top, Bottom
    LC_Nodes = {}
    NX = (EX+1)
    NY = (EY+1)
    A = LX*LY
    for i in range(NX*NY):
            LC_Nodes.update({i : [q*A/(NX*NY), 0, 0]})
    return(LC_Nodes)

# if __name__=='__main__':
    # (N, DGoF) = CCCC(3,3)
    # print(N)
    # print(DGoF)
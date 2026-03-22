import numpy as np

def rectangularmesh(Lx,Ly,ElementsX,ElementsY, inix, iniy):
    #rectangularmesh(float, float, int, int)
    Nx = ElementsX
    Ny = ElementsY
    xx = np.linspace(inix,Lx,Nx+1)
    yy = np.linspace(iniy,Ly,Ny+1)
    nodeCoor = np.zeros(((Nx+1)*(Ny+1),2))
    elementNodes = np.zeros((Nx*Ny,4))
    for i in range(0,len(xx)):
        for j in range(0,len(yy)):
            call=(Ny+1)*(i)+j
            nodeCoor[call,0] = xx[i]
            nodeCoor[call,1] = yy[j]
    d = 0
    for i in range(0,Nx):
        for j in range(0,Ny):
            # elementNodes[d,1] = np.floor(d/Ny)*(Ny+1) + d%Ny
            # elementNodes[d,0] = np.floor(d/Ny)*(Ny+1) + d%Ny + 1
            elementNodes[d,0] = int(np.floor(d/Ny)*(Ny+1) + d%Ny)
            elementNodes[d,1] = int(np.floor(d/Ny)*(Ny+1) + d%Ny + 1)
            elementNodes[d,2] = int(np.floor(d/Ny)*(Ny+1) + d%Ny + Ny + 1)
            elementNodes[d,3] = int(np.floor(d/Ny)*(Ny+1) + d%Ny + Ny + 2)
            d += 1
    return nodeCoor, elementNodes

def RectToNodes(rect, EY):
    NodeTopLeft = (rect%EY) + np.floor(rect/EY)*(EY + 1)
    NodeBottomLeft = (rect%EY) + np.floor(rect/EY)*(EY + 1) + 1
    NodeTopRight = (rect%EY) + (np.floor(rect/EY) + 1)*(EY + 1)
    NodeBottomRight = (rect%EY) + (np.floor(rect/EY) + 1)*(EY + 1) + 1
    NodeTopLeft = int(NodeTopLeft)
    NodeBottomLeft = int(NodeBottomLeft)
    NodeTopRight = int(NodeTopRight)
    NodeBottomRight = int(NodeBottomRight)
    return [NodeTopLeft, NodeBottomLeft, NodeTopRight, NodeBottomRight]

def NodeToDF(Node):
    DFX = 2*Node
    DFY  = 2*Node + 1
    return [DFX, DFY]

if __name__=='__main__':
    nodeCoor , elementNodes = rectangularmesh(60,60,6,6,0,0)
    print('coord')
    print(nodeCoor)
    print('elementnodes')
    print(elementNodes)

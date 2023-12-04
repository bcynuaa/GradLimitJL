from Boundary import Boundary
import numpy as np

if __name__  == '__main__':

    ##* generate inner and outer nodes
    t = np.linspace(0,2*np.pi,20)
    radius = 0.05
    innerNodes = np.c_[radius*np.cos(t)+0.5,radius*np.sin(t)+0.5]
    outerNodes = np.array([[0,0],[1,0],[1,1],[0,1]])

    ##* set spacing(S) and weights(Psi) of inner and outer nodes
    innerSpacing = np.ones(len(innerNodes))*0.01 # set all inner spacing to 0.01
    outerSpacing = np.ones(len(outerNodes))      # ser all outer spacing to 1

    innerNodalPsi = np.ones(len(innerNodes))*0.75     # set all inner nodal psi to 0.5
    outerNodalPsi = np.ones(len(outerNodes))*0.25     # ser all outer linear psi to 0.25

    innerLinearPsi = np.ones(len(innerNodes))*0.75    # same as nodal psi for simplification
    outerLinearPsi = np.ones(len(outerNodes))*0.25    

    ##* Instantiation and initialization
    innerBD = Boundary()
    outerBD = Boundary()

    # set num of inner/outer nodes 
    innerBD.allocate(len(innerNodes))
    outerBD.allocate(len(outerNodes))

    # set position of inner/outer nodes
    innerBD.setNodesPosition(innerNodes)
    outerBD.setNodesPosition(outerNodes)

    # set spacing of inner/outer nodes
    innerBD.setNodalSpace(innerSpacing)
    outerBD.setNodalSpace(outerSpacing)

    # set psi of inner/outer nodes
    innerBD.setNodalPsi(innerNodalPsi)
    innerBD.setLinearPsi(innerLinearPsi)
    outerBD.setNodalPsi(outerNodalPsi)
    outerBD.setLinearPsi(outerLinearPsi)

    # cartesian grid
    X,Y = np.meshgrid(np.linspace(0,1,512),np.linspace(0,1,512))
    XY = np.stack([X,Y],axis=2).reshape(-1,2)

    sumInnerPsiJ = []
    sumInnerPsiI = []
    for xy in XY:
        sumInnerPsiJ.append(innerBD.calSumPsiJ(xy[0],xy[1]))
        sumInnerPsiI.append(innerBD.calSumPsiI(xy[0],xy[1]))

    sumOuterPsiJ = []
    sumOuterPsiI = []
    for xy in XY:
        sumOuterPsiJ.append(outerBD.calSumPsiJ(xy[0],xy[1]))
        sumOuterPsiI.append(outerBD.calSumPsiI(xy[0],xy[1]))

pass



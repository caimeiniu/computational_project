#!/usr/bin/python
import numpy as np
import math
import random
import scipy.spatial
#
# Create perfect lattice with axes [110], [112], [111] for use 
# in constructing a columnar nanocyrstal
#
latticeConstant=3.639087    
nX=2 
nY=40
nZ=72
atomPosition=np.zeros([nX*nY*nZ*2,3],np.float64)
lX=latticeConstant/np.sqrt(2.0)
lY=latticeConstant*np.sqrt(2.0)*np.sin(60.0*np.pi/180.0)
lZ=latticeConstant/np.sqrt(3.0)
stackingVectorX=3.0*lX/2.0/3.0
stackingVectorY=lY/2.0/3.0
#
# FCC lattice is created by an ABC stacking of 111 planes 
#
iSite=0
iStacking=0
for iZ in range(nZ):
  posZ=iZ*lZ
  for iX in range(nX):
    for iY in range(nY):

      posX=lX*iX+iStacking*stackingVectorX
      posY=lY*iY+iStacking*stackingVectorY
      if (posX>=nX*lX) : posX=posX-nX*lX
      if (posY>=nY*lY) : posY=posY-nY*lY
      iSite=iSite+1
      atomPosition[iSite-1,0]=posX
      atomPosition[iSite-1,1]=posY
      atomPosition[iSite-1,2]=posZ

      posX=lX*iX+iStacking*stackingVectorX+lX/2.0
      posY=lY*iY+iStacking*stackingVectorY+lY/2.0
      if (posX>=nX*lX) : posX=posX-nX*lX
      if (posY>=nY*lY) : posY=posY-nY*lY
      iSite=iSite+1
      atomPosition[iSite-1,0]=posX
      atomPosition[iSite-1,1]=posY
      atomPosition[iSite-1,2]=posZ

  iStacking=iStacking+1
  if (iStacking==3) : iStacking=0
#
# Shift lattice so that the center of mass is at the origin
#
numberOfAtoms=iSite
comX=0.0
comY=0.0
comZ=0.0
for iSite in range(numberOfAtoms):
  comX+=atomPosition[iSite,0]
  comY+=atomPosition[iSite,1]
  comZ+=atomPosition[iSite,2]
comX/=numberOfAtoms
comY/=numberOfAtoms
comZ/=numberOfAtoms
for iSite in range(numberOfAtoms):
  atomPosition[iSite,0]-=comX
  atomPosition[iSite,1]-=comY
  atomPosition[iSite,2]-=comZ
#
# Construct nanocrystal consisting of four such grains, first defining
# grain centres and grain orientations
#
sideLengthX=float(nX)*lX
sideLengthY=200.0
sideLengthZ=200.0*math.sqrt(3.0/4.0)
numberOfGrains=4
grainPosition=np.zeros([numberOfGrains,3],np.float64)
grainPosition[0,]=[0.0,0.0,0.0]
grainPosition[1,]=[0.0,sideLengthY/2.0,0.0]
grainPosition[2,]=[0.0,sideLengthY/4.0,sideLengthZ/2.0]
grainPosition[3,]=[0.0,3.0*sideLengthY/4.0,sideLengthZ/2.0]
grainAngle=np.zeros([numberOfGrains],np.float64)
for iGrain in range(numberOfGrains) :
  grainAngle[iGrain]=random.random()*np.pi
#
# Print grain positions
#
print (grainPosition)
#
# Print grain orientation angles
#
print (grainAngle)
#
# Construct NC configuration
#
nanocrystalPosition=[]
nanocrystalGrain=[]
for iGrain in range(numberOfGrains) :
  for iSite in range(numberOfAtoms) :
      posX=grainPosition[iGrain,0]+atomPosition[iSite,0]
      posY=grainPosition[iGrain,1]+atomPosition[iSite,1]*math.cos(grainAngle[iGrain])-atomPosition[iSite,2]*math.sin(grainAngle[iGrain])
      if (posY>=sideLengthY) :  posY-=sideLengthY
      if (posY<0.0) : posY+=sideLengthY
      posZ=grainPosition[iGrain,2]+atomPosition[iSite,1]*math.sin(grainAngle[iGrain])+atomPosition[iSite,2]*math.cos(grainAngle[iGrain])
      if (posZ>=sideLengthZ) :  posZ-=sideLengthZ
      if (posZ<0.0) : posZ+=sideLengthZ
      minDis2=1e6
      for jGrain in range(numberOfGrains) :
          disY=posY-grainPosition[jGrain,1]
          if (disY>=sideLengthY/2.0) : disY-=sideLengthY
          if (disY<-sideLengthY/2.0) : disY+=sideLengthY
          disZ=posZ-grainPosition[jGrain,2]
          if (disZ>=sideLengthZ/2.0) : disZ-=sideLengthZ
          if (disZ<-sideLengthZ/2.0) : disZ+=sideLengthZ
          if (disY**2+disZ**2<minDis2) :
              minDis2=disY**2+disZ**2
              kGrain=jGrain
      if (kGrain==iGrain) :
          nanocrystalPosition.append([posX,posY,posZ])
          nanocrystalGrain.append(iGrain)
#          
# Check for close atoms (remove 1)
#
numberofatoms=len(nanocrystalPosition)
dist=scipy.spatial.distance.pdist(nanocrystalPosition, metric='euclidean')        
distmatrix=scipy.spatial.distance.squareform(dist)
booleanmtx=(distmatrix)<latticeConstant/(2*2**0.5)
problem=[]
for i in range(numberofatoms):
    for j in range(i+1,numberofatoms):
        if (distmatrix[i][j]<latticeConstant/(2*2**0.5)): 
            problem.append(i)
            break
nanocrystalPosition = [i for i in nanocrystalPosition if nanocrystalPosition.index(i) not in problem]
#  
# Write data to lammps file
#
f=open("mdyn_nc.lammps","w+")
f.write("FCC crystal\r\n\r\n")
f.write("%d atoms\r\n\r\n" % (len(nanocrystalPosition)))
f.write("%d atom types\r\n\r\n" % (2))
f.write("%f %f xlo xhi\r\n" % (-sideLengthX/2.0,sideLengthX/2.0))
f.write("%f %f ylo yhi\r\n" % (0.0,sideLengthY))
f.write("%f %f zlo zhi\r\n\r\n" % (0.0,sideLengthZ))
f.write("Atoms\r\n\r\n")
for iSite in range(len(nanocrystalPosition)):
  f.write("%d %d %f %f %f\r\n" % (iSite+1,1,nanocrystalPosition[iSite][0],nanocrystalPosition[iSite][1],nanocrystalPosition[iSite][2]))
f.close()


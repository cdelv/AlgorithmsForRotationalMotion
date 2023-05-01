# encoding: utf-8

#---------------------------------------------------------------------
# Simple (rigid) clump bouncing inside a box
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Carlos Andrés del Valle <cdelv@unal.edu.co> (Vasileios did almost everything =)
# Institute for Multiscale Simulation

# Units: m, N, kg, Pa, sec

#---------------------------------------------------------------------
# Import python libraries and yade modules
from __future__ import print_function
from yade import pack,ymport,export,geom,plot
from builtins import range
import random
import numpy as np
import os
import shutil

#---------------------------------------------------------------------
# Material
matID=O.materials.append(FrictMat(young=1.0e9, poisson=0.3, density=2500.0, frictionAngle=0.0, label='spheresMat'))

# Facet box
id1 = O.bodies.append(geom.facetBox((Vector3(0.5,0.5,0.5)), (Vector3(1.0, 1.0, 1.0)),material=O.materials['spheresMat'], color=(1,0,0), wire=True))

# Create a sphere cloud
sp = pack.SpherePack()
sp.makeClumpCloud((0, 0, 0), (1, 1, 1), 
    [
    pack.SpherePack(
            [
            ((-0.04, -0.04, -0.04), 0.02*sqrt(3)), 
            ((  0.0,   0.0,   0.0), 0.02*sqrt(3)), 
            (( 0.04,  0.04,  0.04), 0.02*sqrt(3))
            ]
        )
    ], 
    periodic=False, num=30, seed=11)

# add the sphere pack to the simulation
sp.toSimulation(material='spheresMat')

# Ensure that all particles have the same kinetic energy
Vel = 6
random.seed(11)
for b in O.bodies:
    if isinstance(b.shape,Clump):
        theta = random.random()
        phi = random.random()
        vv = Vector3(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))
        b.state.vel=Vel*vv
        b.state.angVel=0.5*vv

print("x,y,z,vx,vy,vz,wx,wy,wz,q1,q2,q3,qw")
for b in O.bodies:
    if isinstance(b.shape,Clump):
        v = Vector3(1,0,0)
        x = b.state.pos
        v = b.state.vel
        w = b.state.angVel
        q = b.state.ori
        I = b.state.inertia
        m = b.state.mass
        #print(x[0],",", x[1],",", x[2],",",v[0],",",v[1],",",v[2],",",w[0],",",w[1],",",w[2],",",q[3],",",q[0],",",q[1],",",q[2])
        #print(I[0],",", I[1],",", I[2], ",", m)

#exit()

#---------------------------------------------------------------------
Type = {'Omelyan':False, 'YADE':True}
integrationType='YADE'

RungeKutta = False

# Engines
O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(label='ig2')],
            [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.0, betas=0.0, label='ip2')],
            [Law2_ScGeom_MindlinPhys_Mindlin(label='law')]
        ),
        NewtonIntegrator(damping=0.0, gravity=[0,0,0], label='newton', kinSplit=True,
            exactAsphericalRot=True, Omelyan98=False, Allen89=Type[integrationType]),
        VTKRecorder(fileName='temp-', recorders=['all'], iterPeriod=1, label='vtkRec',dead=True)
]


# Set time step
O.dt=0.7*PWaveTimeStep()
#O.dt=0.1*RayleighWaveTimeStep()
print('dt = ', O.dt)

#---------------------------------------------------------------------
# Plot energies
O.trackEnergy = True
O.step()
plot.plots={'t':('total','kinTrans','kinRot')}

def addPlotData(): 
    plot.addData(t=O.time , total=O.energy.total() , **O.energy )
    plot.saveDataTxt('Test3_' + integrationType + '.txt') # Update results 
    #plot.plot(noShow=True).savefig('Test1_' + integrationType + '.png') # Update plot 

O.engines+=[PyRunner(command='addPlotData()', iterPeriod=200000)]

#------------------------------------------------------------------------------
# Function: Start the user interface if YADE is ran with GUI (i.e. if "yadedaily -n clumpBouncingOnFacetBox.py" is not used)
try:
    from yade import qt
    rndr = yade.qt.Renderer()
    v=qt.View()

except: pass

O.saveTmp()
O.stopAtTime = 250
#O.run()
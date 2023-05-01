from yade import qt
from yade import plot
import csv
import yade.timing
O.timingEnabled=True

# Initial angular velocity
wx0 = 1.0
wy0 = 1.0
wz0 = 1.0

# Inertia tensor components
Ix = 2.0
Iy = 1.0
Iz = 1.0

# Torque
Tx = 0.0

tmax = 20.0
dt = 0.0001
n_samples = 100
n_steps = int(tmax/dt)
vis_steps = int(n_steps/n_samples)

# File data
file = "Yade0.csv"
with open(file, 'w') as f:
    f.write('t,q0,q1,q2,q3,wx,wy,wz\n')

# Create a sphere
Body = utils.sphere([0,0,0], 1.0)

# Make it non spherical and put it to spin
Body.aspherical = True
Body.state.inertia = Vector3(Ix, Iy, Iz)
Body.state.angVel = Vector3(wx0, wy0, wz0)
Body.state.angMom = Vector3(Ix*wx0, Iy*wy0, Iz*wz0)

# Add it to the simulation
O.bodies.append(Body)

# Set time step
O.dt=dt
print('dt = ', O.dt)

# Save simulation
#O.saveTmp()
O.stopAtTime = tmax

# Create the engine
O.engines=[
    ForceResetter(),
    InsertionSortCollider(),
    InteractionLoop(),
    NewtonIntegrator(damping=0.0, gravity=[0,0,0], exactAsphericalRot=True, Omelyan98=False, Allen89=True),
    PyRunner(iterPeriod=vis_steps, command='addplot()')
]

def addplot():
    b = O.bodies[0]
    q = b.state.ori
    w = q.conjugate().toRotationMatrix()*b.state.angVel
    plot.addData(
            wx=w[0],
            wy=w[1],
            wz=w[2],
            t=O.time,
        )
    
    qq = q.conjugate()
    with open(file, 'a') as f:
        writer = csv.writer(f)
        writer.writerow([O.time, qq[0], qq[1], qq[2], qq[3], w[0], w[1], w[2]])


plot.plots = {'t': ('wx', 'wy', 'wz')}
plot.plot(subPlots=False)


# Run the simulation
O.run()

O.wait()

yade.timing.stats()
import numpy as np
import matplotlib.pyplot as plt 

from eqMotionSolver import eqMotion

m0 = 68000
mf = 9714
hturn = 130
d = 5
CD = 0.5
gamma0 = 89.85

tf = 260
step = 10000
term = True #Terminate integration on burnout (True) or coast after burnout (False) 

sol = eqMotion(m0,mf,hturn,d,CD,gamma0,tf,step,term)

print(sol)

#Graphs

## Mass vs time
mvst = plt.figure(num="Mass [kg] vs time [s]")
mvstp = mvst.add_subplot()
mvstp.plot(sol.t,sol.y[4])
mvstp.grid()
mvst.supxlabel("Time [s]")
mvst.supylabel("Mass [kg]")
mvst.suptitle("Mass [kg] vs time [s]")


## Speed vs time
vvst = plt.figure(num="Speed [km/s] vs time [s]")
vvstp = vvst.add_subplot()
vvstp.plot(sol.t,sol.y[0]/1000)
vvstp.grid()
vvst.supxlabel("Time [s]")
vvst.supylabel("Speed [km/s]")
vvst.suptitle("Speed [km/s] vs time [s]")

## Altitude vs downrange
xvsh = plt.figure(num="Altitude [km] vs downrange [km]")
xvshp = xvsh.add_subplot()
xvshp.plot(sol.y[2]/1000,sol.y[3]/1000)
xvshp.grid()
xvsh.supxlabel("Downrange [km]")
xvsh.supylabel("Altitude [km]")
xvsh.suptitle("Altitude [km] vs downrange [km]")

## Flight path angle [deg] vs time [s]
phivst = plt.figure(num="Flight path angle [deg] vs time [s]")
phivstp = phivst.add_subplot()
phivstp.plot(sol.t,np.rad2deg(sol.y[1]))
phivstp.grid()
phivst.supxlabel("Time [s]")
phivst.supylabel("Flight path angle [deg]")
phivst.suptitle("Flight path angle [deg] vs time [s]")


plt.show()
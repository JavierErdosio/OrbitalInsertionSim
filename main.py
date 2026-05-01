import numpy as np
import matplotlib.pyplot as plt 

from eqMotionSolver import eqMotion

hturn = 130
d = 5
CD = 0.5
phi0 = 89.85

tf = 260
step = 10000
term = True #Terminate integration on burnout (True) or coast after burnout (False) 


stages = {"stage1":{ 
            "m0": 68000,
            "mf": 9714,
            "Thrust": 933912,
            "ISP": 390
            }
          }

tComplete = np.array([0])
vComplete = np.array([0])
phiComplete = np.array([np.deg2rad(phi0)])
xComplete = np.array([0])
hComplete = np.array([0])
massComplete = np.array([stages["stage1"]["m0"]])

for i in stages:
    sol = eqMotion(stages[i]["m0"],stages[i]["mf"],stages[i]["Thrust"],stages[i]["ISP"],hturn,d,CD,vComplete[-1],phiComplete[-1],xComplete[-1],hComplete[-1],tf,step,term)
    
    tComplete=np.concatenate((tComplete,sol.t))
    vComplete=np.concatenate((vComplete,sol.y[0]))
    phiComplete=np.concatenate((phiComplete,sol.y[1]))
    xComplete=np.concatenate((xComplete,sol.y[2]))
    hComplete=np.concatenate((hComplete,sol.y[3]))
    massComplete=np.concatenate((massComplete,sol.y[4]))



#Graphs

## Mass vs time
mvst = plt.figure(num="Mass [kg] vs time [s]")
mvstp = mvst.add_subplot()
mvstp.plot(tComplete,massComplete)
mvstp.grid()
mvst.supxlabel("Time [s]")
mvst.supylabel("Mass [kg]")
mvst.suptitle("Mass [kg] vs time [s]")


## Speed vs time
vvst = plt.figure(num="Speed [km/s] vs time [s]")
vvstp = vvst.add_subplot()
vvstp.plot(tComplete,vComplete/1000)
vvstp.grid()
vvst.supxlabel("Time [s]")
vvst.supylabel("Speed [km/s]")
vvst.suptitle("Speed [km/s] vs time [s]")

## Altitude vs downrange
xvsh = plt.figure(num="Altitude [km] vs downrange [km]")
xvshp = xvsh.add_subplot()
xvshp.plot(xComplete/1000,hComplete/1000)
xvshp.grid()
xvsh.supxlabel("Downrange [km]")
xvsh.supylabel("Altitude [km]")
xvsh.suptitle("Altitude [km] vs downrange [km]")

## Flight path angle [deg] vs time [s]
phivst = plt.figure(num="Flight path angle [deg] vs time [s]")
phivstp = phivst.add_subplot()
phivstp.plot(tComplete,np.rad2deg(phiComplete))
phivstp.grid()
phivst.supxlabel("Time [s]")
phivst.supylabel("Flight path angle [deg]")
phivst.suptitle("Flight path angle [deg] vs time [s]")


plt.show()

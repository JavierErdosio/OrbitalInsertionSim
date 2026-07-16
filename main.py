import numpy as np
import matplotlib.pyplot as plt 

from eqMotionSolver import eqMotion
from SEZtoECEF import SEZtoECEF
from orbitalParams import orbitalParams


################# DATA ######################
### Launch Site
Lat = 28   # [deg] (phi)
Lon = -80   # [deg] (theta)
H = 0       # [km] Altitude to mean sea level

### Launch parameters
beta = 50      #[deg] Launch azimuth measured from north and clockwise (90 = east)

hturn = 3800     # Altitude to start performing active controlled variation of flight path angle
hg = 8200        # Altitude to start performing gravity turn
Adotphi = -0.025 # Active controlled variation of flight path angle
d = 3.7        # Rocket diameter
CD = 0.3       # Drag coefficient
phi0 = 89.5    #Initial flight path angle

tf = 500       #Maximum simulation time per stage
step = 10000   #Number of steps per stage
term = True    #Terminate integration on burnout (True) or coast after burnout (False) 

PL = 28*800  # Payload mass

stages = {"stage1":{ 
            "m0": PL + (25600+395700) +  (3900+92670), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (25600+220000) +  (3900+92670),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 7800e3,
            "ISP": 300
            },
          "stage1.1":{ #MaxQ
            "m0": PL + (25600+220000) +  (3900+92670), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (25600+200000) +  (3900+92670),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 7800e3*0.6,
            "ISP": 300
            },
          "stage1.2":{ 
            "m0": PL + (25600+200000) +  (3900+92670), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (25600+25700) +  (3900+92670),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 7800e3,
            "ISP": 300
            },
          "stage2.1":{ #Fairing deployment
            "m0": PL + (0+0) +  (3900+92670),          #Payload + (          0            +            0          ) + (Second stage structural + Second stage propellant)
            "mf": PL + (0+0) +  (3900+92000),              #Payload + (          0            +            0          ) + (Second stage structural +             0          )
            "Thrust": 981e3*.1,
            "ISP": 348
            },
          "stage2.2":{ 
            "m0": PL + (0+0) +  (2900+92000),          #Payload + (          0            +            0          ) + (Second stage structural + Second stage propellant)
            "mf": PL + (0+0) +  (2900+1200),              #Payload + (          0            +            0          ) + (Second stage structural +             0          )
            "Thrust": 981e3*0.94,
            "ISP": 348
            }
          }


############### SOLUTION ####################
tComplete = np.array([0])
vComplete = np.array([0])
phiComplete = np.array([np.deg2rad(phi0)])
xComplete = np.array([0])
hComplete = np.array([0])
massComplete = np.array([stages["stage1"]["m0"]])

tEvents = {}
yEvents = {}

for i in stages:
    sol = eqMotion(stages[i]["m0"],stages[i]["mf"],stages[i]["Thrust"],stages[i]["ISP"],hturn,d,CD,vComplete[-1],phiComplete[-1],xComplete[-1],hComplete[-1],tf,step,term,hg,Adotphi)
    
    tEvents[i] = sol.t_events+tComplete[-1]
    yEvents[i] = sol.y_events

    tComplete=np.concatenate((tComplete,tComplete[-1]+sol.t))
    vComplete=np.concatenate((vComplete,sol.y[0]))
    phiComplete=np.concatenate((phiComplete,sol.y[1]))
    xComplete=np.concatenate((xComplete,sol.y[2]))
    hComplete=np.concatenate((hComplete,sol.y[3]))
    massComplete=np.concatenate((massComplete,sol.y[4]))

#Downrange decomposition
rSouth = -xComplete*np.cos(np.deg2rad(beta))
rEast = xComplete*np.sin(np.deg2rad(beta))

#Speed vector
sinPhi = np.sin(phiComplete)
cosPhi = np.cos(phiComplete)

verticalSpeed = np.multiply(vComplete,sinPhi)

horizontalSpeed = np.multiply(vComplete,cosPhi)

vSouth = -horizontalSpeed*np.cos(np.deg2rad(beta))
vEast = horizontalSpeed*np.sin(np.deg2rad(beta))


#SEZ to ECEF
rSEZ = np.column_stack((rSouth,rEast,hComplete))
vSEZ = np.column_stack((vSouth,vEast,verticalSpeed))

rECEF,vECEF = SEZtoECEF(Lat,Lon,H,rSEZ,vSEZ)

#Orbital parameters at final burnout
r,v,Rp,Ra,h,inc,omega,RAAN,theta,e = orbitalParams(rECEF[-1]/1000,vECEF[-1]/1000)
print(
    "Distance to center:    %.3f [km] \n" \
    "Current speed:         %.3f [km/s] \n" \
    "Periapsis radius:      %.3f [km] \n" \
    "Apoapsis radius:       %.3f [km] \n" \
    "Angular momentum:      %.3f [km^2/s] \n" \
    "Inclination:           %.3f [deg] \n" \
    "Argument of periapsis: %.3f [deg] \n" \
    "RAAN:                  %.3f [deg] \n" \
    "True anomaly:          %.3f [deg] \n" \
    "Eccentricity:          %.3f [-]" %(r,v,Rp,Ra,h,inc,omega,RAAN,theta,e)
)

################# Graphs ######################
## True
tTrue = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 145, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 515]

hTrue = [0, 0.1, 0.5, 1.4, 2.9, 5, 7.7, 11, 15, 20, 24.5, 30, 36.8, 44, 52.6, 58.9, 71, 79.6, 87.5, 95, 102, 108, 114, 120, 125, 129, 133, 137, 140, 143, 145, 147, 149, 150, 152, 152, 153, 153, 153, 153,152,152,151,150,149,149,148,147,146,145,144,143,143,143]

vTrue = [0.000, 0.028, 0.075, 0.133, 0.198, 0.275, 0.342, 0.439, 0.556, 0.708, 0.881, 1.083, 1.333, 1.611, 1.917, 2.167, 2.150, 2.192, 2.235, 2.287, 2.345, 2.406, 2.472, 2.543, 2.616, 2.694, 2.778, 2.867, 2.960, 3.058, 3.160, 3.268, 3.380, 3.498, 3.621, 3.751, 3.887, 4.028, 4.178, 4.335, 4.501, 4.677, 4.863, 5.059, 5.268, 5.489, 5.721, 5.971, 6.235, 6.513, 6.814, 7.181, 7.513, 7.581]

maxQ = [71,0.439]

## Mass vs time
mvst = plt.figure(num="Mass [kg] vs time [s]")
mvstp = mvst.add_subplot()
mvstp.plot(tComplete,massComplete)
for i in tEvents:
    mvstp.plot([tEvents[i][0][0],tEvents[i][0][0]],[0,yEvents[i][0][0][4]],linestyle="dashed",label=i+" Burnout: %.2f"%yEvents[i][0][0][4])
mvstp.grid()
mvstp.legend()
mvst.supxlabel("Time [s]")
mvst.supylabel("Mass [kg]")
mvst.suptitle("Mass [kg] vs time [s]")


## Speed vs time
vvst = plt.figure(num="Speed [km/s] vs time [s]")
vvstp = vvst.add_subplot()
vvstp.plot(tComplete,vComplete/1000)
for i in tEvents:
    vvstp.plot([tEvents[i][0][0],tEvents[i][0][0]],[0,yEvents[i][0][0][0]/1000],linestyle="dashed",label=i+" Burnout: %.2f"%(yEvents[i][0][0][0]/1000))
vvstp.scatter(tTrue,vTrue,label="Falcon 9 16-10-2025")
vvstp.plot([maxQ[0],maxQ[0]],[0,maxQ[1]],linestyle="dashed",label="maxQ Falcon 9 16-10-2025")
vvstp.grid()
vvstp.legend()
vvst.supxlabel("Time [s]")
vvst.supylabel("Speed [km/s]")
vvst.suptitle("Speed [km/s] vs time [s]")

## Altitude vs downrange
xvsh = plt.figure(num="Altitude [km] vs downrange [km]")
xvshp = xvsh.add_subplot()
xvshp.plot(xComplete/1000,hComplete/1000)
for i in tEvents:
    xvshp.plot([yEvents[i][0][0][2]/1000,yEvents[i][0][0][2]/1000],[0,yEvents[i][0][0][3]/1000],linestyle="dashed",label=i+" Burnout. h = %.2f" %(yEvents[i][0][0][3]/1000))
xvshp.grid()
xvshp.legend()
xvsh.supxlabel("Downrange [km]")
xvsh.supylabel("Altitude [km]")
xvsh.suptitle("Altitude [km] vs downrange [km]")


## Altitude vs time
tvsh = plt.figure(num="Altitude [km] vs time [s]")
tvshp = tvsh.add_subplot()
tvshp.plot(tComplete,hComplete/1000)
for i in tEvents:
     tvshp.plot([tEvents[i][0][0],tEvents[i][0][0]],[0,yEvents[i][0][0][3]/1000],linestyle="dashed",label=i+" Burnout: %.2f"%(yEvents[i][0][0][3]/1000))
tvshp.scatter(tTrue,hTrue,label="Falcon 9 16-10-2025")
tvshp.grid()
tvshp.legend()
tvsh.supxlabel("Time [s]")
tvsh.supylabel("Altitude [km]")
tvsh.suptitle("Altitude [km] vs time [s]")

## Flight path angle [deg] vs time [s]
phivst = plt.figure(num="Flight path angle [deg] vs time [s]")
phivstp = phivst.add_subplot()
phivstp.plot(tComplete,np.rad2deg(phiComplete))
for i in tEvents:
    phivstp.plot([tEvents[i][0][0],tEvents[i][0][0]],[0,np.rad2deg(yEvents[i][0][0][1])],linestyle="dashed",label=i+" Burnout: %.2f" %np.rad2deg(yEvents[i][0][0][1]))
phivstp.grid()
phivstp.legend()
phivst.supxlabel("Time [s]")
phivst.supylabel("Flight path angle [deg]")
phivst.suptitle("Flight path angle [deg] vs time [s]")


plt.show()

import numpy as np
import matplotlib.pyplot as plt 

from eqMotionSolver import eqMotion
from TOPtoECEF import SEZtoECEF,ENZtoECEF
from orbitalParams import orbitalParams
from TwoBodySolver import SatPoints


################# DATA ######################
### Launch Site
Lat = -38.96   # [deg] (phi)
Lon = -61.71   # [deg] (theta)
H = 0       # [km] Altitude to mean sea level

### Launch parameters
beta = 185      #[deg] Launch azimuth measured from north and clockwise (90 = east)

hturn = 200    # Altitude to start performing gravity turn
hg = 2200
Adotphi = -0.007 # Active controlled variation of flight path angle
d = 2.5        # Rocket diameter
CD = 0.3       # Drag coefficient
phi0 = 88.66    #Initial flight path angle

tf = 1000       #Maximum simulation time per stage
step = 10000   #Number of steps per stage
term = True    #Terminate integration on burnout (True) or coast after burnout (False) 

PL = 500  # Payload mass

test = 3000

stages = {"stage1":{ 
            "m0": PL + (5600+60000) +  (900+5000), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (5600+20000) +  (900+5000),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 107.25e3*9.81,
            "ISP": 262
            },
            "stage1.1":{ 
            "m0": PL + (5600+20000) +  (900+5000), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (5600+0) +  (900+5000),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 107.25e3*9.81*0.8,
            "ISP": 262
            },
            "stage2.1":{ 
            "m0": PL + (0+0) +  (700+5000), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (0+0) +  (700+4500),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 2975*9.81*0.95,
            "ISP": 317
            },
            "stage2.2":{ 
            "m0": PL + (0+0) +  (700+4500), #Payload + (First stage structural + First stage propellant) + (Second stage structural + Second stage propellant)
            "mf": PL + (0+0) +  (700+test),      #Payload + (First stage structural +            0          ) + (Second stage structural + Second stage propellant)
            "Thrust": 2975*9.81,
            "ISP": 317
            }
          }


############### SOLUTION ####################
tComplete = np.array([0])
vComplete = np.array([0])
phiComplete = np.array([np.deg2rad(phi0)])
xComplete = np.array([0])
hComplete = np.array([H])
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

rECEF,vECEF,rLaunchSite,vLaunchSite = SEZtoECEF(Lat,Lon,H,rSEZ,vSEZ)

#Propagation
rogvog = np.concatenate((rECEF[-1],vECEF[-1]))
time,pos,vel = SatPoints(1000,step,rogvog,PL+700+test,test,2975*9.81*.88,317,"Angular",14)

rECEFpos = np.concatenate((rECEF,pos))
vECEFpos= np.concatenate((vECEF,vel))
timepos = np.concatenate((tComplete,time+tComplete[-1]))


#Orbital parameters at final burnout
r,v,Rp,Ra,h,inc,omega,RAAN,theta,e = orbitalParams(rECEFpos[-1],vECEFpos[-1])
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

altitudeECEF = np.linalg.norm(rECEFpos,axis=1)-rLaunchSite
velECEF = np.linalg.norm(vECEFpos,axis=1)

################# Graphs ######################
## True
tTrue = [0,156,156.7,158,160,720]

hTrue = [0,108,109,112.6,116.9,600]

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
vvstp.plot(tComplete,vComplete/1000,label="SEZ")
vvstp.plot(timepos,velECEF,label="ECI")
for i in tEvents:
    vvstp.plot([tEvents[i][0][0],tEvents[i][0][0]],[0,yEvents[i][0][0][0]/1000],linestyle="dashed",label=i+" Burnout: %.2f"%(yEvents[i][0][0][0]/1000))
vvstp.grid()
vvstp.legend()
vvst.supxlabel("Time [s]")
vvst.supylabel("Speed [km/s]")
vvst.suptitle("Speed [km/s] vs time [s]")

## Altitude vs downrange
xvsh = plt.figure(num="Altitude [km] vs downrange [km]")
xvshp = xvsh.add_subplot()
xvshp.plot(xComplete/1000,hComplete/1000,label="SEZ")
#xvshp.plot(xComplete/1000,altitudeECEF,label="ECEF")
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
tvshp.plot(tComplete,hComplete/1000,label="SEZ")
tvshp.plot(timepos,altitudeECEF,label="ECI")
for i in tEvents:
     tvshp.plot([tEvents[i][0][0],tEvents[i][0][0]],[0,yEvents[i][0][0][3]/1000],linestyle="dashed",label=i+" Burnout: %.2f"%(yEvents[i][0][0][3]/1000))
tvshp.scatter(tTrue,hTrue,label="TII-250 Ejemplo")
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

import numpy as np
from scipy.optimize import fsolve


############### Constants ###############
muEarth = 398600 #[km^3/s^2]
Re = 6378 #Earth radius
hours = 23.93446944 #Sidereal day
seg = hours*3600
omegaEarth = 2*np.pi/seg #Earth angular velocity
g0 = 9.81 #[m/s^2] Earth gravity at sea level

#Assumptions
DeltaVL = 1 #[km/s] DeltaV due to losses
RBurnout = 600 + Re #[km] Rocket position radius measured from Earth center
phi = np.deg2rad(0) #[deg] Flight path angle - 0 due to circular orbit or apogee o perigee
beta = np.deg2rad(90+0) #[deg] Launch azimuth (measured from north and clockwise)

#Launch site
Lat = np.deg2rad(-41) #[deg] Latitude of launch site
Lon = np.deg2rad(-63) #[deg] Longitude of launch site
RL = Re #Earth radius at launch site

#Rocket constants
PL = 18800 #Payload mass

Isp1 = 297.5 #[s] Falcon 9 first stage average
m01 = PL + (25600+395700) +  (3900+92670) #[kg] Initial mass - (Payload + (Dry+Propellant)_s1 + (Dry+Propellant)_s2)
mf1 = PL + (25600+0) +  (3900+92670) #[kg] Final mass (mass at burnout)

Isp2 = 315.5 #[s] Falcon 9 second
m02 = PL + (0+0) +  (3900+92670) #[kg] Initial mass 
mf2 = PL + (0+0) +  (3900+0)   #[kg] Final mass (mass at burnout)


############### Equations ###############
DeltaVD = Isp1*g0*np.log(m01/mf1)+Isp2*g0*np.log(m02/mf2) #[m/s] Designed DeltaV 

DeltaVN = DeltaVD/1000 - DeltaVL #[km/s] DeltaV available to get to orbit (DeltaV needed)

#DeltaVN is a vectorial composition of the DeltaV needed to reach certain altitude plus the velocity at burnout and minus the velocity at launch site 
#DeltaVN = DeltaVPE + DeltaVBurnout - DeltaVLaunch_site

DeltaVPE = ((2*muEarth*(RBurnout-RL))/(RL*RBurnout))**(0.5) #DeltaV due to change in altitude

#Launch due east
DeltaVLaunch_site = omegaEarth*(Re*np.cos(Lat)) #DeltaV due to Earth rotation

vecBurnout = np.array([-np.cos(phi)*np.cos(beta),np.cos(phi)*np.sin(beta),np.sin(phi)]) #Velocity vector direction at burnout

func =  lambda DeltaVBurnout: ((DeltaVBurnout*vecBurnout[0])**2+(DeltaVBurnout*vecBurnout[1]-DeltaVLaunch_site)**2+(DeltaVPE+DeltaVBurnout*vecBurnout[2])**2) - DeltaVN**2

vBurnout = fsolve(func,5000)[0] #Find velocity at burnout


############### Orbital parameters ###############

epsilon = vBurnout**2/2 - muEarth/RBurnout #[km^2/s^2] Specific mechanical energy

a = -muEarth/(2*epsilon) #[km] semimajor axis

Ra = 2*a-RBurnout #[km] Apogee radius

e = (Ra-RBurnout)/(Ra+RBurnout) #[-] Eccentricity

if Lat < 0:
    i = -np.acos(np.sin(beta)*np.cos(Lat)) #[rad] inclination
else:
    i = np.acos(np.sin(beta)*np.cos(Lat)) #[rad] inclination

if e>=0:
    print("Orbit has been reached with the following parameters: \n Burnout velocity = %.3f [km/s] \n Eccentricity = %.3f [-] \n Perigee radius = %.3f [km] \n Apogee radius = %.3f [km] \n Inclination = %.3f [deg]" %(vBurnout,e,RBurnout,Ra,np.rad2deg(i)))

else:
    print("Orbit hasn't been reached (e = %.3f), try reducing payload mass" %(e))
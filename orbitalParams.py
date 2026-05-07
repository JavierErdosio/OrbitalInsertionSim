import numpy as np


muEarth = 398600 #[km^3/s^2]

def orbitalParams(r0,v0):  #[km] and [km/s]

    #1)
    r = np.linalg.norm(r0) #Distance to planet center

    #2)
    v = np.linalg.norm(v0) #Speed 

    #3)
    vr = (v0 @ r0)/r #Radial speed

    #4)
    h_vec = np.cross(r0,v0) #Specific angular momentum vector

    #5)
    h = np.linalg.norm(h_vec) #Specific angular momentum

    #6)
    i = np.acos(h_vec[2]/h) #Inclination

    #7)
    N_vec = np.cross(np.array([0,0,1]),h_vec)

    #8)
    N = np.linalg.norm(N_vec)

    #9)
    if N_vec[1] < 0:
        RAAN = 2*np.pi-np.acos(N_vec[0]/N) #Right ascension of the ascending node
    else:
        RAAN = np.acos(N_vec[0]/N)

    #10)
    e_vec = 1/muEarth*((v**2 - muEarth/r)*r0 - r*vr*v0) #Eccentricity vector

    #11)
    e = np.linalg.norm(e_vec) #Eccentricity

    #12)
    if e_vec[2]<0:
        omega = 2*np.pi-np.acos((N_vec @ e_vec)/(N*e)) #Argument of perigee
    else:
        omega = np.acos((N_vec @ e_vec)/(N*e))

    #13)
    theta = np.acos((e_vec @ r0)/(e*r)) #True anomaly

    #14)
    a = r*(1+e*np.cos(theta))/(1-e**2) #Semimajor axis

    Rp = a*(1-e) #Periapsis
    Ra = a*(1+e) #Apoapsis

    return r,v,Rp,Ra,h,np.rad2deg(i),np.rad2deg(omega),np.rad2deg(RAAN),np.rad2deg(theta),e
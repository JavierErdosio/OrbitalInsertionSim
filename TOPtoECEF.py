import numpy as np



def SEZtoECEF(Lat,Lon,H,rSEZ,vSEZ):
    #Earth
    Re = 6378*1000   #[m] Earth Radius
    omegaEarth = 2*np.pi/(23.93446944*3600) #[rad/s]

    #Topocentric
    f = 0.003353 # Earth Oblateness

    Lat = np.deg2rad(Lat)
    Lon = np.deg2rad(Lon)

    Rphi = Re/(1-(2*f-f**2)*np.sin(Lat)**2)**0.5

    Rc = Rphi + H
    Rs = ((1-f)**2)*Rphi + H

    x0 = Rc*np.cos(Lat)
    z0 = Rs*np.sin(Lat)

    LaunchSiteSpeedSEZ= [0,omegaEarth*x0,0]

    #ECEF
    X = x0*np.cos(Lon)
    Y = x0*np.sin(Lon)
    Z = z0

    LaunchSitePos = np.array([X,Y,Z])

    ########### SEZ to ECEF ##############
    sin_phi = np.sin(Lat)
    cos_phi = np.cos(Lat)
    sin_theta = np.sin(Lon)
    cos_theta = np.cos(Lon)

    RotSEZtoECEF = np.array([
        [ sin_phi*cos_theta, -sin_theta, cos_phi*cos_theta],
        [ sin_phi*sin_theta,  cos_theta, cos_phi*sin_theta],
        [-cos_phi,          0.0,     sin_phi]
        ])

    LaunchSiteSpeed = np.matvec(RotSEZtoECEF, LaunchSiteSpeedSEZ)

    rECEF = LaunchSitePos + np.matvec(RotSEZtoECEF, rSEZ)
    vECEF = LaunchSiteSpeed + np.matvec(RotSEZtoECEF, vSEZ)

    return rECEF/1000, vECEF/1000,np.linalg.norm(LaunchSitePos)/1000,LaunchSiteSpeed/1000


def ENZtoECEF(Lat,Lon,H,rSEZ,vSEZ):
    #Earth
    Re = 6378*1000   #[m] Earth Radius
    omegaEarth = 2*np.pi/(23.93446944*3600) #[rad/s]

    #Topocentric
    f = 0.003353 # Earth Oblateness

    Lat = np.deg2rad(Lat)
    Lon = np.deg2rad(Lon)

    Rphi = Re/(1-(2*f-f**2)*np.sin(Lat)**2)**0.5

    Rc = Rphi + H
    Rs = ((1-f)**2)*Rphi + H

    x0 = Rc*np.cos(Lat)
    z0 = Rs*np.sin(Lat)

    LaunchSiteSpeedENZ= [omegaEarth*x0,0,0]

    #ECEF
    X = x0*np.cos(Lon)
    Y = x0*np.sin(Lon)
    Z = z0

    LaunchSitePos = np.array([X,Y,Z])

    ########### SEZ to ECEF ##############
    sin_phi = np.sin(Lat)
    cos_phi = np.cos(Lat)
    sin_theta = np.sin(Lon)
    cos_theta = np.cos(Lon)

    RotENZtoECEF = np.array([
        [-sin_theta, -sin_phi * cos_theta,  cos_phi * cos_theta],
        [ cos_theta, -sin_phi * sin_theta,  cos_phi * sin_theta],
        [ 0.0,         cos_phi,               sin_phi             ]
    ])

    LaunchSiteSpeed = np.matvec(RotENZtoECEF, LaunchSiteSpeedENZ)

    rECEF = LaunchSitePos + np.matvec(RotENZtoECEF, rSEZ)
    vECEF = LaunchSiteSpeed + np.matvec(RotENZtoECEF, vSEZ)

    return rECEF/1000, vECEF/1000,np.linalg.norm(LaunchSitePos)/1000,LaunchSiteSpeed/1000
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 

Re = 6378*1000 #[m] Earth radius
g0 = 9.81 #[m/s^2]
rho0 = 1.225 #kg/m^3
m0 = 68000 #[kg] initial mass

def func(t,y):

    v,phi,x,h,m = y[0],y[1],y[2],y[3],y[4]    

    g = g0/(1+h/Re)**2 #[m/s^2] Gravity

    mburnout = 9714 #[kg] Mass at burnout 

    if m <= mburnout: 
        T= 0
        dotm = 0
    
    else:
        T = 1.4 * m0 * g0 #[N] Thrust
        Isp = 390 #[s]
        dotm = -T/(Isp*g0) #Mass flow rate

    rho=rho0*np.exp(-h/7500) #[kg/m^3] Density
    
    d = 5 #[m]
    A = np.pi*d**2/4 #[m^2] Frontal area

    CD = 0.5 #[-] Drag coefficient

    D = 0.5*rho*v**2*A*CD #[N] Drag force

    hturn = 130 #[m] Altitude to start performing gravity turn

    if h <= hturn:
        dotphi = 0 #[rad/s] Variation of flight path angle
        dotv = T/m-D/m-g #[m/s^2] Acceleration
        dotx = 0 #[m/s] Horizontal speed
        doth = v #[m/s] Vertical speed
    
    else:
        if v < 1e-3:
            dotphi = 0 #[rad/s] Variation of flight path angle
        else:
            dotphi = -(1/v) * (g-v**2/(Re+h))*np.cos(phi) #[rad/s] Variation of flight path angle

        dotv = T/m-D/m-g*np.sin(phi) #[m/s^2] Acceleration
        dotx = (Re/(Re+h))*v*np.cos(phi) #[m/s] Horizontal speed
        doth = v*np.sin(phi) #[m/s] Vertical speed
 

    return dotv,dotphi,dotx,doth,dotm



sol = solve_ivp(func,[0,260],
                [0,np.deg2rad(89.85),0,0,m0]
                ,method='DOP853' #Explicit Runge-Kutta method of order 8.
                ,rtol=1e-8
                ,atol=1e-10
                ,t_eval=np.arange(0,260,step=260/10000)
                )

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
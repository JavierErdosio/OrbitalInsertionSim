import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 

Re = 6378*1000 #[m] Earth radius
g0 = 9.81 #[m/s^2]
rho0 = 1.225 #kg/m^3
m0 = 68000 #[kg]

def func(t,y):

    v,phi,x,h,m = y[0],y[1],y[2],y[3],y[4]    

    g = g0/(1+h/Re)**2 #[]

    if m <= 9714:
        T= 0
        dotm = 0
    
    else:
        T = 1.4 * m0 * g0
        Isp = 390
        dotm = -T/(Isp*g0)

    rho=rho0*np.exp(-h/7500) #[]
    
    d = 5 #[m]
    A = np.pi*d**2/4 #[m^2]

    CD = 0.5 #[-]

    D = 0.5*rho*v**2*A*CD

    hturn = 130

    if h <= hturn:
        dotphi = 0
        dotv = T/m-D/m-g
        dotx = 0
        doth = v
    
    else:
        if v < 1e-3:
            dotphi = 0
        else:
            dotphi = -(1/v) * (g-v**2/(Re+h))*np.cos(phi)

        dotv = T/m-D/m-g*np.sin(phi)
        dotx = (Re/(Re+h))*v*np.cos(phi)
        doth = v*np.sin(phi)
 

    return dotv,dotphi,dotx,doth,dotm


#print(func(0,[0,np.deg2rad(90),0,0,68000]))
sol = solve_ivp(func,[0,260],
                [0,np.deg2rad(89.85),0,0,m0]
                ,method='DOP853'
                ,rtol=1e-8
                ,atol=1e-10
                ,t_eval=np.arange(0,260,step=260/10000)
                )

mvst = plt.figure(num="Mass [kg] vs time [s]")
mvstp = mvst.add_subplot()
mvstp.plot(sol.t,sol.y[4])
mvstp.grid()
mvst.supxlabel("Time [s]")
mvst.supylabel("Mass [kg]")

vvst = plt.figure(num="Speed [km/s] vs time [s]")
vvstp = vvst.add_subplot()
vvstp.plot(sol.t,sol.y[0]/1000)
vvstp.grid()
vvst.supxlabel("Time [s]")
vvst.supylabel("Speed [km/s]")

xvsh = plt.figure(num="Downrange [km] vs altitude [km]")
xvshp = xvsh.add_subplot()
xvshp.plot(sol.y[2]/1000,sol.y[3]/1000)
xvshp.grid()
xvsh.supxlabel("Downrange [km]")
xvsh.supylabel("Altitude [km]")

phivst = plt.figure(num="Flight path angle [deg] vs time [s]")
phivstp = phivst.add_subplot()
phivstp.plot(sol.t,np.rad2deg(sol.y[1]))
phivstp.grid()
phivst.supxlabel("Time [s]")
phivst.supylabel("Flight path angle [deg]")


plt.show()
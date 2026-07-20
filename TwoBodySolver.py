import numpy as np
from numpy import cos,sin
from scipy.integrate import odeint
from scipy.integrate import solve_ivp


#Datos generales
muTierra = 398600 #[km^3/s^2]
R_e = 6378 #[km]
f = 0.003353
g0 = 9.81



def SatPoints(seconds,steps,rogvog,mt,mp,Thrust,Isp):
    tf = seconds #Tiempo [s] 
    pasos = steps
    mburnout = mt-mp

    #Resolucion
    y0 = np.concatenate((rogvog,[mt])) #Vector de estado
    t = np.linspace(0, tf, pasos)

    def ec(t,y):
        X,Y,Z,Vx,Vy,Vz,m = y
        
        pos = np.linalg.norm([X,Y,Z])
        vel= np.linalg.norm([Vx,Vy,Vz])

        if m <= mburnout: 
            T= 0
            dotm = 0
        
        else:
            T = Thrust #[N] Thrust
            dotm = -T/(Isp*g0) #Mass flow rate

        # Ecuaciones de estado
        Ax = -muTierra*(X)/pos**3+Vx/vel*T/m/1000
        Ay = -muTierra*(Y)/pos**3+Vy/vel*T/m/1000
        Az = -muTierra*(Z)/pos**3+Vz/vel*T/m/1000

        return([Vx,Vy,Vz,Ax,Ay,Az,dotm])

    sol = solve_ivp(ec,
                    [0,tf],
                    y0,t_eval=t,
                    method='DOP853',
                    rtol=1e-10,
                    atol=1e-14)

    
    x = sol.y[0]
    y = sol.y[1]
    z = sol.y[2]
    vx = sol.y[3]
    vy = sol.y[4]
    vz = sol.y[5]

    points = np.column_stack((x,y,z))
    vel= np.column_stack((vx,vy,vz))

    return t,points,vel

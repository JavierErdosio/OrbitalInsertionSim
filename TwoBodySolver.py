import numpy as np
from numpy import cos,sin
from scipy.integrate import solve_ivp


#Datos generales
muTierra = 398600 #[km^3/s^2]
R_e = 6378 #[km]
f = 0.003353
g0 = 9.81



def SatPoints(seconds,steps,rogvog,mt,mp,Thrust,Isp,Modes,ang=0,OOF=True):
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

        hvec = np.cross(np.array([X,Y,Z]),np.array([Vx,Vy,Vz]))
        h = np.linalg.norm(hvec)

        verNorm = hvec/h
        verRad = np.array([X,Y,Z])/pos
        
        vecTan = np.cross(verNorm,verRad)
        verTan = vecTan/np.linalg.norm(vecTan)

        velVers = np.array([Vx/vel,Vy/vel,Vz/vel])


        if m <= mburnout: 
            T= 0
            dotm = 0
                
        else:
            T = Thrust #[N] Thrust
            dotm = -T/(Isp*g0) #Mass flow rate

        if abs(pos-6378-580) <20:
            mode = "Velocity"
        else:
            mode = Modes


        if mode == "Velocity":
            Vers = velVers
        
        else:
            if mode == "Horizontal":
                Vers = verTan
            
            elif mode == "Angular":
                theta = np.radians(ang)
                Vers = cos(theta) * verTan + sin(theta) * verRad

        

        # Ecuaciones de estado
        Ax = -muTierra*(X)/pos**3+Vers[0]*T/m/1000
        Ay = -muTierra*(Y)/pos**3+Vers[1]*T/m/1000
        Az = -muTierra*(Z)/pos**3+Vers[2]*T/m/1000


        return([Vx,Vy,Vz,Ax,Ay,Az,dotm])

    def outOfProp(t,y):
            return y[6]-mburnout
        
    outOfProp.terminal = OOF

    sol = solve_ivp(ec,
                    [0,tf],
                    y0,t_eval=t,
                    method='DOP853',
                    rtol=1e-10,
                    atol=1e-14,
                    events=outOfProp)

    
    x = sol.y[0]
    y = sol.y[1]
    z = sol.y[2]
    vx = sol.y[3]
    vy = sol.y[4]
    vz = sol.y[5]

    points = np.column_stack((x,y,z))
    vel= np.column_stack((vx,vy,vz))

    return sol.t,points,vel

import numpy as np
from scipy.integrate import solve_ivp

Re = 6378*1000 #[m] Earth radius
g0 = 9.81 #[m/s^2]
rho0 = 1.225 #kg/m^3


def eqMotion(m0,mburnout,hturn,d,CD,gamma0,tf,steps,term):  #Initial mass [kg], Mass at burnout [kg], Altitude to start performing gravity turn [m], rocket diameter [m], Drag coefficient [-], Initial flight path angle [deg], Final simultation tife [] 
    def func(t,y):

        v,phi,x,h,m = y[0],y[1],y[2],y[3],y[4]    

        g = g0/(1+h/Re)**2 #[m/s^2] Gravity

        if m <= mburnout: 
            T= 0
            dotm = 0
        
        else:
            T = 1.4 * m0 * g0 #[N] Thrust
            Isp = 390 #[s]
            dotm = -T/(Isp*g0) #Mass flow rate

        rho=rho0*np.exp(-h/7500) #[kg/m^3] Density
        
        A = np.pi*d**2/4 #[m^2] Frontal area

        D = 0.5*rho*v**2*A*CD #[N] Drag force

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

    def outOfProp(t,y):
        return y[4]-mburnout
    
    outOfProp.terminal = term


    sol = solve_ivp(func,[0,tf],
                    [0,np.deg2rad(gamma0),0,0,m0]
                    ,method='DOP853' #Explicit Runge-Kutta method of order 8.
                    ,rtol=1e-8
                    ,atol=1e-10
                    ,t_eval=np.arange(0,tf,step=tf/steps)
                    ,events=outOfProp
                    )
    
    return sol


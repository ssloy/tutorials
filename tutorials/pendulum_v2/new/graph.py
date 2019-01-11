#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

#[time, refcurrent, current, x, theta,target] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)
[time, refcurrent, current, x, theta] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum run")
ax1.set_xlabel('Time, sec')

#ax1.plot(time, refcurrent, color='green',   label='reference current, A')
#ax1.plot(time, current,    color='red',    label='measured current, A')
#ax1.plot(time, x,          color='blue',    label='cart position, m')
ax1.plot(time, theta,      color='magenta', label='pendulum angle, rad')


hatX = x[0]
hatQ = theta[0]
hatLPX = 0.
hatLPQ = 0.
hatSigmaX = 0
hatSigmaQ = 0

hatDX = 0.
hatDQ = 0.


hatXm = x[0]
hatQm = theta[0]
hatLPXm = 0.
hatLPQm = 0.

hatDXm = 0.
hatDQm = 0.


x_dot = [0]
theta_dot = [0]
x_dot_observer_m = [0]
q_dot_observer_m = [0]

x_dot_observer = [0]
q_dot_observer = [0]


for i in range(len(time)-1):
    dt = (time[i+1]-time[i])

    x_dot.append((x[i+1]-x[i])/dt)
    theta_dot.append((theta[i+1]-theta[i])/dt)

#             current_ref = (int16_t)((f + (v==0L ? (f>0L?4690L:-3350L) : (v>0L?4690L:-3350L)))/13L);

    u = current[i+1]*12.3392
    if abs(x_dot_observer_m[i])>.001:
        if (x_dot_observer_m[i]>0):
            u = u - 4.3323
        else:
            u = u + 3.0295


    mesQ = theta[i+1]
    mesX = x[i+1]
    cosQ = np.cos(mesQ)
    sinQ = np.sin(mesQ)


# PEBO
    diffXm = 0.83150999999998020939528942108154*hatLPXm - 50.0*hatXm + 50.0*mesX - (97.103000000002793967723846435547*hatLPQm*cosQ)/np.sqrt(68340.0 - 13637.0*cosQ*cosQ);
    diffQm = 90.0*mesQ - 90.0*hatQm + (1114.6000000000931322574615478516*hatLPQm) / np.sqrt(68340.0 - 13637.0*cosQ*cosQ);
    diffLPXm = 332.60999999998603016138076782227*mesX - 332.60999999998603016138076782227*hatXm + 0.83150999999998020939528942108154*u + (0.000000000000000000000000002117582368135664566385455895796*cosQ* (18343113252400251770071909662720.0*hatQm - 18343113252400251770071909662720.0*mesQ))/np.sqrt(68340.0 - 13637.0*cosQ*cosQ);
    diffLPQm = -(0.000000000000000000000000000033087224502119758849772748371813* (5.0530186655360592957378995945472e33*hatQm - 5.0530186655360592957378995945472e33*mesQ - 41639056032387337505159667253248.0*sinQ + 2934770911164556041706136403968.0*u*cosQ))/np.sqrt(68340.0 - 13637.0*cosQ*cosQ);
    hatDXm = 0.83150999999998020939528942108154*hatLPXm - (9710300.0*hatLPQm*cosQ)/np.sqrt(683400000000000.0 - 136370000000000.0*cosQ*cosQ);
    hatDQm = (1114.6000000000931322574615478516*hatLPQm)/np.sqrt(68340.0 - 13637.0*cosQ*cosQ);

    hatXm   = hatXm   + dt*diffXm
    hatQm   = hatQm   + dt*diffQm
    hatLPXm = hatLPXm + dt*diffLPXm
    hatLPQm = hatLPQm + dt*diffLPQm

    x_dot_observer_m.append(hatDXm)
    q_dot_observer_m.append(hatDQm)


# HGO6
    A = 1./(cosQ*cosQ - 4.803)
    diffX = hatDX + 384.615*(mesX - hatX);  
    diffQ = hatDQ + 384.615*(mesQ - hatQ);  
    diffDX = hatSigmaX + 29585.799*(mesX - hatX) - (4.178*u + .449*hatDQ*hatDQ*sinQ - 9.810*cosQ*sinQ  ) * A; 
    diffDQ = hatSigmaQ + 29585.799*(mesQ - hatQ) + (9.295*u*cosQ - 104.828*sinQ + hatDQ*hatDQ*cosQ*sinQ) * A; 
    diffSigmaX = 1820664.543*(mesX - hatX); 
    diffSigmaQ = 1820664.543*(mesQ - hatQ); 

#   diffX = hatDX + 2500.*(mesX - hatX);  
#   diffQ = hatDQ + 2500.*(mesQ - hatQ);  
#   diffDX = hatSigmaX + 1250000.*(mesX - hatX) - (4.178*u + .449*hatDQ*hatDQ*sinQ - 9.810*cosQ*sinQ  ) * A; 
#   diffDQ = hatSigmaQ + 1250000.*(mesQ - hatQ) + (9.295*u*cosQ - 104.828*sinQ + hatDQ*hatDQ*cosQ*sinQ) * A; 
#   diffSigmaX = 500000000.0*(mesX - hatX); 
#   diffSigmaQ = 500000000.0*(mesQ - hatQ); 

    hatX  = hatX  + dt*diffX;
    hatQ  = hatQ  + dt*diffQ;
    hatDX = hatDX + dt*diffDX;
    hatDQ = hatDQ + dt*diffDQ;
    hatSigmaX = hatSigmaX + dt*diffSigmaX;
    hatSigmaQ = hatSigmaQ + dt*diffSigmaQ;



    x_dot_observer.append(hatDX)
    q_dot_observer.append(hatDQ)

#ax1.plot(time, x_dot_observer         , color='red',   label='xdot HGO6')
#ax1.plot(time, x_dot_observer_m, color='pink', label='xdot PEBO')

#ax1.plot(time, q_dot_observer, color='green', label='thetadot HGO6')
#ax1.plot(time, q_dot_observer_m, color='black', label='thetadot PEBO')

ax1.plot(time, refcurrent         , color='blue',   label='refcurrent')
ax1.plot(time, current         , color='red',   label='current')

ax1.legend()
plt.show()


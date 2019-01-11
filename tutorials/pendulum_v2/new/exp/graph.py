#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

[time, refcurrent, current, x, theta, target] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum run")
ax1.set_xlabel('Time, sec')

#ax1.plot(time, refcurrent, color='green',   label='reference current, A')
#ax1.plot(time, current,    color='red',    label='measured current, A')
ax1.plot(time, x,          color='blue',    label='cart position, m')
ax1.plot(time, theta,      color='magenta', label='pendulum angle, rad')
ax1.plot(time, target,     color='red', label='pendulum angle, rad')


hatX = x[0]
hatQ = theta[0]
hatLPX = 0.
hatLPQ = 0.

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

#    diffDXm = 30.8*np.sign(mesX - hatXm) - (5.5784376938378675056094284022284e37*u + 7.0288314942357130570678797868078e36*hatDQm*hatDQm*sinQ - 1.5794273385254500758861270837494e38*cosQ*sinQ)/ (1.610017674337869598252932807084e37*cosQ*cosQ - 8.0681593206050990965838229134771e37);
#    diffDQm = 85.8*np.sign(mesQ - hatQm) + (5.111167220120220946834707324076e34* (100.0*u*cosQ - 1418.8317102217761075610980014972*sinQ + 12.6*hatDQm*hatDQm*cosQ*sinQ))/(6.4400706973514783930117312283358e35*cosQ*cosQ - 3.2272637282420396386335291653908e36);
#    diffXm = hatDXm + 7.9372539331937717715048472609178*np.sign(mesX - hatXm)*np.sqrt(np.abs(hatXm - mesX))
#    diffQm = hatDQm + 13.247641299491770282146064090439*np.sign(mesQ - hatQm)*np.sqrt(np.abs(hatQm - mesQ))
#
#    hatXm   = hatXm   + dt*diffXm
#     hatQm   = hatQm   + dt*diffQm
#    hatDXm = hatDXm + dt*diffDXm
#    hatDQm = hatDQm + dt*diffDQm


    x_dot_observer_m.append(hatDXm)
    q_dot_observer_m.append(hatDQm)




#    A = np.sqrt(683 - 136*cosQ*cosQ)
#    hatDQ = (111*hatLPQ) / A
#    hatDX = 0.8315*hatLPX - (9*hatLPQ*cosQ) / A
#    diffX = hatDX - 50*hatX + 50*mesX
#    diffQ = 90*mesQ - 90*hatQ + hatDQ
#    diffLPX = 332*mesX - 332*hatX + 0.8315*u + (3884*cosQ*(hatQ - mesQ)) / A
#    diffLPQ = (16719*(mesQ - hatQ) + 137*sinQ - 9*u*cosQ) / A

#    hatX   = hatX   + dt*diffX
#    hatQ   = hatQ   + dt*diffQ
#    hatLPX = hatLPX + dt*diffLPX
#    hatLPQ = hatLPQ + dt*diffLPQ

#    A = 1/(cosQ*cosQ - 5.011)
#    diffDX = 30.8*np.sign(mesX - hatX) - (3.465*u + .437*hatDQ*hatDQ*sinQ - 9.810*cosQ*sinQ) * A
#    diffDQ = 85.8*np.sign(mesQ - hatQ) + (7.936*u*cosQ - 112.601*sinQ + hatDQ*hatDQ*cosQ*sinQ) * A
#    diffX = hatDX +  7.937*np.sign(mesX - hatX)*np.sqrt(np.abs(hatX - mesX))
#    diffQ = hatDQ + 13.248*np.sign(mesQ - hatQ)*np.sqrt(np.abs(hatQ - mesQ))

    A = 1./(cosQ*cosQ - 5.011)
    diffDX = 2222.0*mesX - 2222.0*hatX - (3.465*u + .437*hatDQ*hatDQ*sinQ -   9.810*cosQ*sinQ) * A
    diffDQ = 2222.0*mesQ - 2222.0*hatQ + (7.936*u*cosQ - 112.601*sinQ + hatDQ*hatDQ*cosQ*sinQ) * A
    diffX = hatDX - 100.0*hatX + 100.0*mesX
    diffQ = hatDQ - 100.0*hatQ + 100.0*mesQ

    hatX   = hatX   + dt*diffX
    hatQ   = hatQ   + dt*diffQ
    hatDX = hatDX + dt*diffDX
    hatDQ = hatDQ + dt*diffDQ


    x_dot_observer.append(hatDX)
    q_dot_observer.append(hatDQ)

#ax1.plot(time, x_dot_observer         , color='red',   label='xdot')
#ax1.plot(time, x_dot_observer_m, color='pink', label='xdot observer')

#ax1.plot(time, q_dot_observer, color='green', label='thetadot')
#ax1.plot(time, q_dot_observer_m, color='black', label='qdot observer')

ax1.legend()
plt.show()


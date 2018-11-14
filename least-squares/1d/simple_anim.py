import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()

#x = np.arange(0., 1., 0.01)
n = 100
f = [0.] * n 
g = [0.] * n 
lock = [0, n*57//100, n-1]
lval = [1., 2., 1.]
for i in range(len(lock)):
    f[lock[i]] = lval[i]

A = np.matrix(np.zeros((n, n)))
b = np.matrix(np.zeros((n, 1)))
for i in range(n):
    if i in lock:
        b[i, 0] = lval[lock.index(i)]
        A[i, i] = 1.
    elif i<n-1 and i>0:
        A[i, i-1] = -.5
        A[i, i]   =  1.
        A[i, i+1] = -.5

#print(A,b)
#f = ((np.linalg.inv(A)*b).transpose().tolist()[0]) 

lines = [ax.plot(f, drawstyle='steps-mid')[0], ax.plot(g, drawstyle='steps-mid')[0]]

def init():  # only required for blitting to give a clean slate.
    for line in lines:
        line.set_ydata([np.nan] * n)
    return lines

def animate(i):
    global f, g

    for i in range(n):
        if i in lock: continue
        g[i] = (f[i-1]+f[i+1])/2. - f[i]

    f = [sum(x) for x in zip(f, g)]
    lines[0].set_ydata(f)  # update the data.
#    lines[1].set_ydata([1.*v for v in g])  # update the data.
    return lines

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=np.arange(0, 100), interval=1, blit=True, save_count=50)

#ani.save('line.gif', dpi=80, writer='imagemagick')

# To save the animation, use e.g.
#
#ani.save("movie.mp4")
#
# or
#
#from matplotlib.animation import FFMpegWriter
#writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#ani.save("movie.mp4", writer=writer)

plt.show()

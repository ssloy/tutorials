import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()

#x = np.arange(0., 1., 0.01)
n = 8
f = [0.] * n
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

print(A,b)
#f = ((np.linalg.inv(A)*b).transpose().tolist()[0]) 

lines = [ax.plot(range(n), f, drawstyle='steps-mid')[0], ax.text(0.05, 0.05, "Iteration #0", transform=ax.transAxes, fontsize=14,bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))]
plt.draw()
#ax.grid()

def animate(iteration):
    global f, n, lock
    if (0==iteration):
        return lines

    norm = 0
    for i in range(n):
#       if (0==i):
#           f[0] = f[1]
#           continue
#       if (i==n-1):
#           f[n-1] = f[n-2]
#           continue

        if i in lock: continue
        val = (f[i-1]+f[i+1]+11./n**2)/2.
#        val = (f[i-1]+f[i+1])/2.
        norm += np.abs(f[i]-val)
        f[i] = val

    if 1 and (n<1000 and norm<1e-2):
        f = [val for val in f for _ in (0, 1)]
        n *= 2
        lock = [0, n*57//100, n-1]
        for i in range(len(lock)):
            f[lock[i]] = lval[i]

    print(iteration, norm)

    lines[0].set_data(range(n), f)  # update the data.
    lines[1].set_text("Iteration #" + str(iteration))
    plt.draw()
    ax.relim()
    ax.autoscale_view(False,True,False)
    return lines

ani = animation.FuncAnimation(fig, animate, frames=np.arange(0, 150), interval=100, blit=False, save_count=50)
ani.save('line.gif', dpi=80, writer='imagemagick')

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

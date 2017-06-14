#! /usr/local/bin/python2.7
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from particle import particle
import sys,os

dt = 0.01 #s
T = 1000
scale = 1.0e-5 #cm
N = 40
appr_speed_at_20C = 175400 #cm/s


big_R = np.empty(shape=N,dtype=particle)
big_N = np.empty(shape=N,dtype=particle)

for part in range(N):
    q = np.random.choice([-1.0,1.0])
    if q == 1.0:
        big_R[part] = particle(1.6726219e-27,1.0e-13, 1.0, np.random.rand(3)*scale, np.random.rand(3)*appr_speed_at_20C, scale) # protons
        big_N[part] = particle(1.6726219e-27,1.0e-13, 1.0, np.random.rand(3)*scale, np.random.rand(3)*appr_speed_at_20C, scale) 

    else:
        #big_R[part] = particle(9.10938356e-31,0.0, -1.0, np.random.rand(3)*scale, np.random.rand(3)*appr_speed_at_20C, scale) # electrons
        #big_N[part] = particle(9.10938356e-31,0.0, -1.0, np.random.rand(3)*scale, np.random.rand(3)*appr_speed_at_20C, scale)
        big_R[part] = particle(1.6726219e-27,1.0e-13, -1.0, np.random.rand(3)*scale, np.random.rand(3)*appr_speed_at_20C, scale) # positrons*
        big_N[part] = particle(1.6726219e-27,1.0e-13, -1.0, np.random.rand(3)*scale, np.random.rand(3)*appr_speed_at_20C, scale)

def kinetic():
    return sum([.5*part.m()*((part.v()*part.v()).sum()) for part in big_R])

for t in range(T):
    fn1 = np.zeros((N,3))
    for part in range(N):
        for i in range(N):
            if i != part:
                fn1[part,:] += big_R[part].force(big_R[i])
        big_N[part].set_r(big_R[part].r() + big_R[part].v()*dt + 0.5*fn1[part,:]*dt*dt)
        
    for part in range(N):
        fn2 = np.zeros(3)
        for i in range(N):
            if i != part:
                fn2 += big_N[part].force(big_N[i])
        #if big_R[part].r().any() > scale:
        #    fn2 = -fn2
        big_N[part].set_v(big_R[part].v() + 0.5*(fn1[part,:]+fn2)*dt)
    big_R = big_N

np.save("./big_R_"+str(N)+"N_"+str(T)+"t_EM.npy",big_R)

#big_R = np.load("./big_R_2000N_500t_GR.npy")
#exit(0)

X = [part.r()[0] for part in big_R]
Y = [part.r()[1] for part in big_R]
Z = [part.r()[2] for part in big_R]
fig = plt.figure()
plt.cool()
'''
colorStar = np.zeros((scale,scale,scale))
for part in big_R:
    colorStar[int(part.r()[0]),int(part.r()[1]),int(part.r()[2])] += 1.0
colorStar = colorStar/colorStar.max()
'''
color = []
for part in big_R:
    #color.append(colorStar[int(part.r()[0]),int(part.r()[1]),int(part.r()[2])])
    color.append((part.q()+1.0)/2.0)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X,Y,Z,c=color)
ax.axis("tight")
#ax.set_xlim3d(0.0,scale)
#ax.set_ylim3d(0.0,scale)
#ax.set_zlim3d(0.0,scale)
plt.title("Gas in a Box, "+str(N)+" Particle T="+str(T)+", Gravity + EM")

plt.savefig(os.path.join('.',"big_R_"+str(N)+"N_"+str(T)+"t_EM.png"),bbox_inches='tight')
plt.close()
#plt.show()


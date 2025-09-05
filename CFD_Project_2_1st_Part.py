# -*- coding: utf-8 -*-
"""
Created on Sun May  1 14:51:06 2022

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib


Mark=["--v","--p","--*","--+","--x"]
vm=1
T=10
v=np.array([vm])
t=np.array([0])
k=-1
leg=[]
L=[]
N=[]
for a in [5,10,20,40,80]:
    dt=1/a
    k=k+1
    vm=1
    v=np.array([vm])
    t=np.array([0])
    n=T/dt
    s=0
    for i in range(1,int(n+1)):
        v1=vm-dt/2*vm
        vmp=-dt*v1+dt*np.sin(i*dt)+vm
        v=np.append(v,[vmp])
        t1=0+i*dt
        t=np.append(t,[t1])
        Ee=np.sin((i-1)*dt)/2+3*np.exp(-(i-1)*dt)/2+np.cos((i-1)*dt)/2-np.sin(i*dt)/2-3*np.exp(-i*dt)/2-np.cos(i*dt)/2
        E=abs(1/dt*Ee-vmp)
        s=s+E
        vm=vmp
    L.append(np.log10(s/n**2))
    N.append(np.log10(n))    
    r=random.randint(8,15)
    a=str(a)
    leg.append(r'$\Delta t= \frac{1}{%s} $'%a)
    ax=plt.gca()
    ax.plot(t,v,Mark[k],markevery=int(len(t)/r))
    ax.set_xlabel('$t$')
    ax.set_ylabel(r'$ \nu (t) $')
    ax.grid(True)

tex=np.linspace(0,10,1000)
vex=3/2*np.exp(-tex)+1/2*np.sin(tex)-1/2*np.cos(tex)
ax.plot(tex,vex,'--o',markevery=int(len(tex)/10))
leg.append("Exact")    
ax.legend(leg)
plt.savefig('./figs/Solution.eps',format='eps')
plt.show()

ax2=plt.gca()
ax2.plot(N,L,'-o')
ax2.legend([r'$ log(Error) $'])
ax2.set_xlabel('$log(N)$')
ax2.set_ylabel(r'$ log(Error) $')
slope=(L[0]-L[-1])/(N[0]-N[-1])
matplotlib.rc('text',usetex=True)
text1='Slope=%g'%slope
ax2.text(N[0], L[0], text1)
ax2.grid(True)
plt.savefig('./figs/L1.eps',format='eps')
plt.show()



# -*- coding: utf-8 -*-
"""
Created on Sun May  1 17:57:25 2022

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import interactive
import matplotlib as mpl
import random

interactive(True)

#%matplotlib qt
def IC(T,l,t):
    import numpy as np
    for i in range(0,T.shape[1]):
        T[0][i]=0
    return T
def BC(i,l,t,T):
    for i in range(1,T.shape[0]):
        T[i][0]=np.sin(2*np.pi*(i)*t/(T.shape[0]-1))
    return T


def SecUPW(T,j,m,u,l,t,v,p):
    import numpy as np
    dx=l/T.shape[1]
    dt=t/(T.shape[0]-1)
    if p==2:
        if j>=2:
            A=-u/(2*dx)*(3*T[m-1][j]-4*T[m-1][j-1]+T[m-1][j-2])
           
        elif j==0:
            A=-u/dx*(T[m-1][j]-(2*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1))-T[m-1][j]))
            #A=-u/(dx*2)*(3*T[m-1][j]-4*(2*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1))-T[m-1][j])+(4*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1)))-2*T[m-1][j]-T[m-1][j])
        else:
           A=-u/dx*(T[m-1][j]-T[m-1][j-1])
           # A=-u/(2*dx)*(3*T[m-1][j]-4*T[m-1][j-1]+(2*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1))-T[m-1][j-1]))          
    elif p==1:
        if j>=1:
            A=-u/dx*(T[m-1][j]-T[m-1][j-1])
        else:
            A=-u/dx*(T[m-1][j]-2*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1))+T[m-1][j])
    elif p==3:
        if j>=1 and j< (T.shape[1]-1):
            #A=-u/(2*dx)*(T[m-1][j+1]-T[m-1][j-1])
            A=-u/(2*dx)*(T[m-1][j+1]-2*T[m-1][j]+T[m-1][j-1])
        elif j==0:
            #A=-u/(2*dx)*(T[m-1][j+1]-2*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1))+T[m-1][j])
            A=-u/(2*dx)*(T[m-1][j+1]-2*T[m-1][j]+2*np.sin(2*np.pi*(m-v*1)*t/(T.shape[0]-1))-T[m-1][j])
        else:
            A=-u/(2*dx)*(3*T[m-1][j]-4*T[m-1][j-1]+T[m-1][j-2])            
    return A 


def RK2(T,m,u,l,t,p):
    import numpy as np
    dx=l/T.shape[1]
    dt=t/(T.shape[0]-1) 
    T1=T
    T3=np.zeros((T.shape[0],T.shape[1]))
    for j in range(0,T.shape[1]):
        a=T[m-1][j]+SecUPW(T,j,m,u,l,t,1,p)*dt/2
        T3[m-1][j]=a
    T4=np.zeros((T.shape[0],T.shape[1]))
    for j in range(0,T.shape[1]):
        b=T[m-1][j]+SecUPW(T3,j,m,u,l,t,1/2,p)*dt
        T4[m][j]=b       
    return T4[m,:]        

def RK4(T,m,u,l,t,p):
    import numpy as np
    dx=l/T.shape[1]
    dt=t/(T.shape[0]-1) 
    T1=T
    T3=np.zeros((T.shape[0],T.shape[1]))
    for j in range(0,T.shape[1]):
        a=T[m-1][j]+SecUPW(T,j,m,u,l,t,1,p)*dt/4
        T3[m-1][j]=a
    T4=np.zeros((T.shape[0],T.shape[1]))
    for j in range(0,T.shape[1]):
        a=T[m-1][j]+SecUPW(T3,j,m,u,l,t,3/4,p)*dt/3
        T4[m-1][j]=a
    T5=np.zeros((T.shape[0],T.shape[1]))
    for j in range(0,T.shape[1]):
        a=T[m-1][j]+SecUPW(T4,j,m,u,l,t,2/3,p)*dt/2
        T5[m-1][j]=a
    T6=np.zeros((T.shape[0],T.shape[1]))
    for j in range(0,T.shape[1]):
        b=T[m-1][j]+SecUPW(T5,j,m,u,l,t,1/2,p)*dt
        T6[m][j]=b       
    return T6[m,:]   


            
nn=2
dt=1/(160*nn)
CFL=0.25
u=2
dx=u*dt/CFL
t=1
l=1
nt=t/dt
nx=1/dx
T=np.zeros((int(nt)+1,int(nx)))
T=IC(T,l,t)
TT=np.linspace(0,1,int(nt)+1)
XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx))
s=0
v=1
w=0


w=0
L1=np.array([])
NN=np.array([])
for nn in [1,2]:
    dt=1/(160*nn)
    CFL=0.25
    u=2
    dx=u*dt/CFL
    t=1
    l=1
    nt=t/dt
    nx=1/dx
    T2=np.zeros((int(nt)+1,int(nx)))
    T2=IC(T2,l,t)
    TT=np.linspace(0,1,int(nt)+1)
    XX2=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx))
    s=0
    v=1
    w=0
    NN=np.append(NN,nx)       
    for k in range(1,int(nt)+1):        
        T2[k,:]=RK2(T2,k,u,l,t,2)  
    for i in range(1,int(nt)+1):
        for j in range(0,int(nx)):
            if XX2[j]<=2*TT[i]:
               # Ee=(1/(dt*dx)*(-np.sin(np.pi*(2*dt*i - dx*j)) + np.sin(((2*i - 2)*dt - j*dx)*np.pi) + np.sin(((-j - 1)*dx + 2*i*dt)*np.pi) - np.sin(np.pi*(2*dt*i - dx*j - 2*dt - dx))))/(2*np.pi**2) 
                Ee=-1/dx*(np.cos(dx*j*np.pi - 2*np.pi*TT[i]) - np.cos(dx*j*np.pi + dx*np.pi - 2*np.pi*TT[i]))/np.pi
        
                E=abs(Ee-T2[i][j])
            else:
                Ee=0
                E=abs(Ee-T2[i][j])
            if i==int(nt):
                w=w+dx*E
    L1=np.append(L1,w)      

ax4=plt.gca()
ax4.plot(np.log10(NN),np.log10(L1))
ax4.grid()
ax4.set_xlabel(r'$Log(N)$')
ax4.set_ylabel(r'$Log(L_1)$')          
Ls=np.log10(L1)
NNs=np.log10(NN)
slope=(Ls[0]-Ls[1])/(NNs[0]-NNs[1])
mpl.rc('text',usetex=True)
text1='Slope=%g'%slope
ax4.text(NNs[0], Ls[0], text1)
plt.savefig('./figs/L1secU.eps',format='eps')




                



XX,TT=np.meshgrid(XX2,TT)
fig, ax=plt.subplots(subplot_kw={"projection":'3d'})

ax.plot_surface(TT,XX2,T2)
plt.show()

        
TT=np.linspace(0,1,1000)
XXe=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),1000) 

Te=np.zeros((1000,1000))

for i in range(0,1000):
    for j in range(0,1000):
        if XXe[j]<=2*TT[i]:
            Te[i][j]=np.sin(2*np.pi*(TT[i]-XXe[j]/2))
        else:
            Te[i][j]=0
            
                                  
XX,TT=np.meshgrid(XXe,TT)
fig, ax=plt.subplots(subplot_kw={"projection":'3d'})

ax.plot_surface(TT,XX,Te,label=r'$Exact Solution$')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel(r'$T_{Exact}$')
fake2Dline = mpl.lines.Line2D([0],[0], linestyle="none", c='b', marker = 'o')
ax.legend([fake2Dline], [r'Exact  Solution'], numpoints = 1)
plt.show()


Mark=["--v","--p","--*","--+","--x"]
leg=[]
kk=-1
for nn in [50,40,35]:
    kk=kk+1
    t=0.25
    dt=t/nn
    nx=40
    CFL=0.25
    u=2
    l=1
    dx=l/40
    nt=t/dt
    nt=nn
    T22=np.zeros((int(nt)+1,int(nx)))
    T22=IC(T22,l,t)
    TT=np.linspace(0,1,int(nt)+1)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx))
    s=0
    v=1
    w=0
    NN=np.append(NN,nx)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx)) 
    a=str(nn)
    leg.append('Time Steps='+a)  
    r=random.randint(8,20)
    for k in range(1,int(nt)+1):        
        T22[k,:]=RK2(T22,k,u,l,t,2)  
    ax3=plt.gca()
    ax3.plot(XX,T22[int(nt),:],Mark[kk],markevery=int(len(XX)/r))
    ax3.grid(True)
    ax3.legend(leg)
    ax3.set_xlabel(r'$x$')
    ax3.set_ylabel(r'$T$')
plt.savefig('./figs/T0.025.eps',format='eps')    
plt.show()



leg=[]
kk=-1
for nx in [40,80]:
    kk=kk+1
    t=1
    dt=t/nn
    CFL=0.25
    u=2
    l=1
    dx=l/nx
    dt=CFL*dx/u
    nt=t/dt
    nt=int(nt)
    T21=np.zeros((int(nt)+1,int(nx)))
    T21=IC(T21,l,t)
    TT=np.linspace(0,1,int(nt)+1)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx))
    NN=np.append(NN,nx)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx)) 
    a=str(nx)
    leg.append('1st order Upwind, Second order Runge-Kutta No. of CVs='+a)  
    r=random.randint(8,20)
    for k in range(1,int(nt)+1):        
        T21[k,:]=RK2(T21,k,u,l,t,1)  
    ax4=plt.gca()
    ax4.plot(XX,T21[int(nt),:],Mark[kk],markevery=int(len(XX)/r))
    
    
leg.append(r'2nd order Upwind, Second order Runge-Kutta No. of CVs=40')
leg.append(r'Exact Solution')
ax4.plot(XX2,T2[-1,:],Mark[kk+1],markevery=int(len(XX2)/r)) 
ax4.plot(XXe,Te[-1,:]) 
ax4.grid(True)
ax4.legend(leg,fontsize=8)
ax4.set_xlabel(r'$x$')
ax4.set_ylabel(r'$T$')
plt.savefig('./figs/T1firstupsecrung.eps',format='eps')    
plt.show()

leg=[]
kk=-1
for CFL in [0.1,0.125,0.5]:
    kk=kk+1
    t=1
    u=2
    l=1
    nx=40
    dt=CFL*dx/u
    dx=l/nx
    nt=t/dt
    nt=int(nt)
    T2c=np.zeros((int(nt)+1,int(nx)))
    T2c=IC(T2c,l,t)
    TT=np.linspace(0,1,int(nt)+1)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx))
    NN=np.append(NN,nx)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx)) 
    a=str(nx)
    #leg.append('1st order Upwind, Second order Runge-Kutta No. of CVs='+a)  
    r=random.randint(8,20)
    for k in range(1,int(nt)+1):        
        T2c[k,:]=RK2(T2c,k,u,l,t,3)  
    a=str(CFL)
    leg.append(r'$CFL=$'+a)
    ax5=plt.gca()
    ax5.grid(True)
    ax5.plot(XX,T2c[int(nt),:],Mark[kk],markevery=int(len(XX)/r))
    ax5.legend(leg)
    ax5.set_xlabel('$x$')
    ax5.set_ylabel('$T$')
plt.savefig('./figs/central.eps',format='eps')    
plt.show()

kk=-1
NN=np.array([])
L14=np.array([])
leg=[]
for nn in [40,80]:
    CFL=0.25
    u=2
    t=1
    l=1
    kk=kk+1
    dt=CFL*dx/u
    nt=t/dt
    nx=nn
    T24=np.zeros((int(nt)+1,int(nx)))
    T24=IC(T24,l,t)
    T222=np.zeros((int(nt)+1,int(nx)))
    T222=IC(T24,l,t)
    TT=np.linspace(0,1,int(nt)+1)
    XX=np.linspace(l/(2*nx),l/(2*nx)*(2*nx-1),int(nx))
    s=0
    v=1
    w=0
    NN=np.append(NN,nx)       
    for k in range(1,int(nt)+1):        
        T24[k,:]=RK4(T24,k,u,l,t,2)
        T222[k,:]=RK2(T24,k,u,l,t,2)
    for i in range(1,int(nt)+1):
        for j in range(0,int(nx)):
            if XX[j]<=2*TT[i]:
               # Ee=(1/(dt*dx)*(-np.sin(np.pi*(2*dt*i - dx*j)) + np.sin(((2*i - 2)*dt - j*dx)*np.pi) + np.sin(((-j - 1)*dx + 2*i*dt)*np.pi) - np.sin(np.pi*(2*dt*i - dx*j - 2*dt - dx))))/(2*np.pi**2) 
                Ee=-1/dx*(np.cos(dx*j*np.pi - 2*np.pi*TT[i]) - np.cos(dx*j*np.pi + dx*np.pi - 2*np.pi*TT[i]))/np.pi
        
                E=abs(Ee-T24[i][j])
            else:
                Ee=0
                E=abs(Ee-T24[i][j])
            if i==int(nt):
                w=w+dx*E
        L14=np.append(L1,w)
    a=str(nx)
    leg.append('RK4 and Second Order Upwind No. of CVs='+a)
    leg.append('RK2 and Second Order Upwind No. of CVs='+a)
    ax7=plt.gca()
    ax7.plot(XX,T24[int(nt),:],Mark[kk],markevery=int(len(XX)/r))
    ax7.plot(XX,T222[int(nt),:],Mark[kk+1],markevery=int(len(XX)/r))
ax7.plot(XXe,Te[-1,:])  
ax7.set_xlabel(r'$x$')
ax7.set_ylabel(r'$T$') 
leg.append('Exact Solution')    
ax7.legend(leg,fontsize=8)
ax7.grid(True) 
plt.savefig('./figs/RK4RK2plots.eps',format='eps')    
plt.show()
ax6=plt.gca()
ax6.plot(np.log10(NN),np.log10(L1))
ax6.grid()
ax6.set_xlabel(r'$Log(N)$')
ax6.set_ylabel(r'$Log(L_1)$')          
Ls=np.log10(L14)
NNs=np.log10(NN)
slope=(Ls[0]-Ls[1])/(NNs[0]-NNs[1])
mpl.rc('text',usetex=True)
text1='Slope=%g'%slope
ax6.text(NNs[0], Ls[0], text1)    
plt.savefig('./figs/L14thrung.eps',format='eps')











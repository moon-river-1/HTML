import numpy as np
import matplotlib.pyplot as plt


xmax=1e-15
dx=1e-17
dt=1e-19
X = np.arange(0,xmax+dx,dx)
M=len(X)-1
timemax=1e-15


alpha=2.76e4 #此处的alpha应该取1还是实际值？
tau_q=0.03e-12
a=alpha/tau_q
tau=1.0

landa=2*dt/(dx*dx*(1-2*tau))

#初始条件
T=np.zeros(M+1)
T[0]=10000
Tt=np.zeros(M+1)
Tt0=np.zeros(M+1)
g1 = T*0.5
g2= T*0.5


#分布函数初始化
g0=Tt+a*landa*T
g1=-a*landa*T/2
g2=-a*landa*T/2


time=0
while time<timemax:
    #源项
    F1=np.zeros(M)
    h1=-dt*(g0+g1+g2)/(3*tau_q)
    
    #碰撞
    g0=(1-1/tau)*g0+(g0+g1+g2+a*landa*T)/tau+h1
    g1=(1-1/tau)*g1-a*landa*T/(2*tau)+h1
    g2=(1-1/tau)*g2-a*landa*T/(2*tau)+h1
    
    #迁移
    for i in range(1,M+1):
        g1[M-i+1] = g1[M-i]
        g2[i-1] = g2[i]
    
    #边界条件
    T[0]=10000
    Tt[0]=0
    g0[0]=Tt[0]+a*landa*T[0]
    g1[0]=-a*landa*T[0]/2
    g2[0]=-a*landa*T[0]/2
    
    
    g0[M]=g0[M-1]
    g1[M]=g1[M-1]
    g2[M]=g2[M-1]
    
    #导数及原函数的计算
    Tt=g0+g1+g2
    T=T+dt*Tt+0.5*(Tt-Tt0)*dt

    time=time+dt
    print(time)
    Tt0=Tt

#图形
fig = plt.figure()
plt.plot(X,T,color='r')
plt.show()
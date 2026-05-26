import numpy as np
from pylab import * #carrega as funções e números básicos (pi, sin, cos, etc..)
from scipy.special import gamma
from matplotlib import pyplot as plt

# import timeit

# start = time.time()

#pm878297623br
# Definição dos parâmetros
ti = 0.0 #tempo inicial
tf = 25 #tempo final
auxT=int(tf-ti)
dt = 0.01
t = np.arange(ti, tf, dt)


n = int(ceil((tf-ti)/dt))
u1 = np.zeros((n),float)
v1 = np.zeros((n),float)
u2 = np.zeros((n),float)
v2 = np.zeros((n),float)


#valores iniciais
u1[0] = 1
v1[0] = 0
u2[0] = 1
v2[0] = 0

#parametros do modelo

m= 1
k = 1
lamb=0.5

d1 = 0.85 # ordem da derivada fracionária.
agamma1 = gamma(d1)


# Resolução do modelo fracionário

for j in range(n-1):
    u1[j+1]=u1[0]+(v1[j])*dt**d1/(d1*agamma1)
    v1[j+1]=v1[0]-((k/m)*u1[j])*dt**d1/(d1*agamma1)
    for i in range(j):
        adt1=((t[j+1]-t[i])**(d1)-(t[j+1]-t[i]-dt)**(d1))/(d1*agamma1)
        u1[j+1] += +(v1[i])*adt1
        v1[j+1] += -((k/m)*u1[i])*adt1
        
# Resolução do modelo interio com atrito
for j in range(n-1):
    u2[j+1]=u2[j]+(v2[j])*dt
    v2[j+1]=v2[j]-((k/m)*u2[j]+(lamb/m)*v2[j])*dt

#Gráficos
plt.figure(figsize=(16,8))
plt.plot(t, u1,'C2', linewidth=3)
plt.plot(t, u2,'C1', linewidth=3)
plt.legend(['$\\alpha=0.85$','$\\alpha=1$'], fontsize=20)
plt.xlabel('$t$', fontsize=20)
plt.ylabel('$x(t)$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20) 
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.grid(True)
plt.savefig('osc.png', format='png')
plt.show()

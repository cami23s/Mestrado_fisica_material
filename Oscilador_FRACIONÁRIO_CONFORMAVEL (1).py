import numpy as np
from pylab import * #carrega as funções e números básicos (pi, sin, cos, etc..)
from scipy.special import gamma
from matplotlib import pyplot as plt

# import timeit

# start = time.time()

#pm878297623br
# Definição dos parâmetros
ti = 0.0 #tempo inicial
tf = 20 #tempo final
auxT=int(tf-ti)
dt = 0.001
t = np.arange(ti, tf, dt)


n = int(ceil((tf-ti)/dt))
u1 = np.zeros((n),float)
v1 = np.zeros((n),float)
u2 = np.zeros((n),float)
v2 = np.zeros((n),float)
u3 = np.zeros((n),float)
v3 = np.zeros((n),float)


#valores iniciais
u1[0] = 1
v1[0] = 0
u2[0] = 1
v2[0] = 0
u3[0] = 1
v3[0] = 0

#parametros do modelo

m= 1
k = 1
lamb=0.0

d1 = 0.5 # ordem da derivada fracionária.
d2 = 0.75 # ordem da derivada fracionária.


# Resolução dos modelos
for j in range(n-1):
    tj=(j+1)*dt
    u1[j+1]=u1[j]+tj**(d1-1)*(v1[j])*dt
    v1[j+1]=v1[j]-(k/m)*u1[j]*tj**(d1-1)*dt
    u2[j+1]=u2[j]+tj**(d2-1)*(v2[j])*dt
    v2[j+1]=v2[j]-(k/m)*u2[j]*tj**(d2-1)*dt
    u3[j+1]=u3[j]+(v3[j])*dt
    v3[j+1]=v3[j]-((k/m)*u3[j]+(lamb/m)*v3[j])*dt
        

print(tj)
print(tj**(d1-1))
#Gráficos
plt.figure(figsize=(16,8))
plt.plot(t, u1,'--C1', linewidth=3)
plt.plot(t, u2,'--C2', linewidth=3)
plt.plot(t, u3, linewidth=3)
plt.legend(['$\\alpha=0.5$', '$\\alpha=0.75$','$\\alpha=1$'], fontsize=20, loc='lower right', framealpha=1)
plt.xlabel('$t$', fontsize=20)
plt.ylabel('$x(t)$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20) 
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.grid(True)
plt.savefig('oscconf.png', format='png')
plt.show()

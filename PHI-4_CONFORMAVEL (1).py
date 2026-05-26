import numpy as np
from pylab import * #carrega as funções e números básicos (pi, sin, cos, etc..)
from scipy.special import gamma
from matplotlib import pyplot as plt

# import timeit

# start = time.time()

#pm878297623br
# Definição dos parâmetros
ti = 0.0 #tempo inicial
tf = 100 #tempo final
auxT=int(tf-ti)
dt = 0.00001
t = np.arange(ti, tf, dt)


n = int(ceil((tf-ti)/dt))
u1 = np.zeros((n),float)
v1 = np.zeros((n),float)
u2 = np.zeros((n),float)
v2 = np.zeros((n),float)
u3 = np.zeros((n),float)
v3 = np.zeros((n),float)
u4 = np.zeros((n),float)
v4 = np.zeros((n),float)
p1 = np.zeros((n),float)
p2 = np.zeros((n),float)
p3 = np.zeros((n),float)
p4 = np.zeros((n),float)
H1 = np.zeros((n),float)
H2 = np.zeros((n),float)
H3 = np.zeros((n),float)
H4 = np.zeros((n),float)


#valores iniciais
u1[0] = 1
v1[0] = 0
u2[0] = 1
v2[0] = 0
u3[0] = 1
v3[0] = 0
u4[0] = 1
v4[0] = 0

#parametros do modelo

m= 1.0
lamb=6.0

d1 = 0.7 # ordem da derivada fracionária.
d2 = 1.3 # ordem da derivada fracionária.


# Resolução dos modelos
p3[0]=v3[0]
p4[0]=v4[0]
H1[0]=p1[0]**2/2+m**2*u1[0]**2/2+lamb*u1[0]**4/24
H2[0]=p2[0]**2/2+m**2*u2[0]**2/2+lamb*u2[0]**4/24
H3[0]=p3[0]**2/2+m**2*u3[0]**2/2+lamb*u3[0]**4/24
H4[0]=p4[0]**2/2+m**2*u4[0]**2/2

for j in range(n-1):
    tj=(j+1)*dt
    u1[j+1]=u1[j]+tj**(d1-1)*(v1[j])*dt
    v1[j+1]=v1[j]-(m**2*u1[j]+(lamb/6.0)*u1[j]**3)*tj**(d1-1)*dt
    u2[j+1]=u2[j]+tj**(d2-1)*(v2[j])*dt
    v2[j+1]=v2[j]-(m**2*u2[j]+(lamb/6.0)*u2[j]**3)*tj**(d2-1)*dt
    u3[j+1]=u3[j]+(v3[j])*dt
    v3[j+1]=v3[j]-(m**2*u3[j]+(lamb/6.0)*(u3[j])**3)*dt
    u4[j+1]=u4[j]+(v4[j])*dt
    v4[j+1]=v4[j]-(m**2*u4[j])*dt
    p1[j+1]=tj**(1-d1)*v1[j+1]
    p2[j+1]=tj**(1-d2)*v2[j+1]
    p3[j+1]=v3[j+1]
    p4[j+1]=v4[j+1]
    H1[j+1]=p1[j+1]**2/2+m**2*u1[j+1]**2/2+lamb*u1[j+1]**4/24
    H2[j+1]=p2[j+1]**2/2+m**2*u2[j+1]**2/2+lamb*u2[j+1]**4/24
    H3[j+1]=p3[j+1]**2/2+m**2*u3[j+1]**2/2+lamb*u3[j+1]**4/24
    H4[j+1]=p4[j+1]**2/2+m**2*u4[j+1]**2/2
    


#Gráficos
plt.figure(figsize=(20,8))
plt.plot(t, u1,'--C1', linewidth=3)
plt.plot(t, u2,'--C2', linewidth=3)
plt.plot(t, u3,'C3', linewidth=3)
plt.plot(t, u4, linewidth=3)
plt.legend(['$\\alpha=0.7$', '$\\alpha=1.3$','$\\alpha=1$','oscilador harmônico'], fontsize=20,loc='lower right', framealpha=1)
plt.xlabel('$t$', fontsize=20)
plt.ylabel('$\\phi(t)$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20) 
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.grid(True)
#plt.savefig('phi4.png', format='png')
plt.show()


#Gráficos
plt.figure(figsize=(20,8))
plt.plot(t, H1,'--C1', linewidth=3)
plt.plot(t, H2,'--C2', linewidth=3)
plt.plot(t, H3,'C3', linewidth=3)
plt.plot(t, H4, linewidth=3)
plt.legend(['$\\alpha=0.7$', '$\\alpha=1.3$','$\\alpha=1$','oscilador harmônico'], fontsize=20,loc='upper left', framealpha=1)
plt.xlabel('$t$', fontsize=20)
plt.ylabel('$H(t)$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20) 
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.grid(True)
#plt.savefig('Hphi4.png', format='png')
plt.show()


fig, ax = plt.subplots(figsize=(10,5.0))

ax.set_aspect('equal')

ax.plot(p1, u1,'--C1', linewidth=3)
ax.plot(p2, u2,'--C2', linewidth=3)
ax.plot(p3, u3,'C3', linewidth=3)
ax.plot(p4, u4, linewidth=3)
ax.set_xlabel('$\\pi$', fontsize=20)
ax.set_ylabel('$\\phi$', fontsize=20)
ax.tick_params(length=10, width=2, labelsize=20)
ax.set_xticks([-4, -3, -2, -1, 0, 1, 2, 3, 4])
ax.set_yticks([-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5])
ax.set_xlim(-4.5, 4.5)
ax.set_ylim(-1.1, 2.6)
ax.grid(True)
ax.legend(['$\\alpha=0.7$', '$\\alpha=1.3$','$\\alpha=1$','oscilador harmônico'], fontsize=14,loc='upper center')
fig.savefig('xpphi4.png', format='png')

plt.show()


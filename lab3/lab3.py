import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
np.seterr(divide='ignore',invalid='ignore')

def potential(x1,y1,x2,y2,m):
	G=6.672e-11
	V=G*(m/np.sqrt((x1-x2)**2+(y1-y2)**2))
	return V
def grav(x1,y1,x2,y2,m):
	G=6.672e-11
	ax=-G*(m*(x1-x2)/np.sqrt((x1-x2)**2+(y1-y2)**2)**3)
	ay=-G*(m*(y1-y2)/np.sqrt((x1-x2)**2+(y1-y2)**2)**3)
	return ax,ay
me=5.9736e24
mm=7.349e22
r_em=384.5e6

x=np.linspace(-4e7,4e8,200)
y=np.linspace(-4e7,2e8,200)

x_a=np.linspace(-4e7,4e8,30)
y_a=np.linspace(-4e7,2e8,30)

x_mesh,y_mesh=np.meshgrid(x,y)
xa_mesh,ya_mesh=np.meshgrid(x_a,y_a)

V1=potential(x_mesh,y_mesh,0,0,me)
V2=potential(x_mesh,y_mesh,r_em*np.cos(10/180*np.pi),r_em*np.sin(10/180*np.pi),mm)
labels=np.floor(np.logspace(0,9))
plt.contour(x_mesh,y_mesh,V1+V2,labels)


ax_e,ay_e=grav(xa_mesh,ya_mesh,0,0,me)
ax_m,ay_m=grav(xa_mesh,ya_mesh,r_em*np.cos(10/180*np.pi),r_em*np.sin(10/180*np.pi),mm)


ax=ax_m+ax_e
ay=ay_m+ay_e
a_norm=np.sqrt(ax**2+ay**2)

plt.quiver(xa_mesh,ya_mesh,ax/a_norm,ay/a_norm,0.2)
plt.show()
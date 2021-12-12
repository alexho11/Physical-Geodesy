import numpy as np
import matplotlib.pyplot as plt

filename = 'Gravimetrie_Hohele.dat'
x,g,h= np.loadtxt(filename,skiprows=1,unpack=True)
# plt.figure()
# plt.subplot(211)
# plt.plot(x,g)
# plt.xlabel('Profilkoordinate x [m]')
# plt.ylabel('g [mgal]')
# plt.subplot(212)
# plt.plot(x,h,color='orange')
# plt.xlabel('Profilkoordinate x [m]')
# plt.ylabel('Hoehe [m]')
# plt.show()

# Korrektur
rho=2700
G=6.674e-11
g_c=g+0.3*h-0.1119*h
# plt.figure()
# plt.plot(x,g_c)
# plt.xlabel('Profilkoordinate x [m]')
# plt.ylabel('g [mgal]')
# plt.title('Schwere nach Korrektur')
# plt.show()
# print(g_c)

# Stoerung
def disturb(x,R,Z,X,rho):
	r=np.sqrt((X-x)**2+(Z+R)**2)
	G=6674e-11
	az=-4/3*np.pi*G*rho*R**3/r**2*np.sin(Z+R/r)
	return az
s1=disturb(600,1890,500,600,2700)
print(s1)

s2=disturb(600,4000,500,600,1000)
print(s2)

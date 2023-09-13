from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Coefficients 	Numerical values
C00 = 1402.388 	
C01 = 5.03830 	
C02 = -5.81090E-2 	
C03 = 3.3432E-4 	
C04 = -1.47797E-6 	
C05 = 3.1419E-9 	
C10 = 0.153563 	
C11 = 6.8999E-4 	
C12 = -8.1829E-6 	
C13 = 1.3632E-7 	
C14 = -6.1260E-10 	
C20 = 3.1260E-5 	
C21 = -1.7111E-6 	
C22 = 2.5986E-8 	
C23 = -2.5353E-10 	
C24 = 1.0415E-12 	
C30 = -9.7729E-9 	
C31 = 3.8513E-10 	
C32 = -2.3654E-12 	
A00 = 1.389 	        
A01 = -1.262E-2
A02 = 7.166E-5   
A03 = 2.008E-6   
A04 = -3.21E-8   
A10 = 9.4742E-5  
A11 = -1.2583E-5 
A12 = -6.4928E-8 
A13 = 1.0515E-8  
A14 = -2.0142E-10
A20 = -3.9064E-7 
A21 = 9.1061E-9  
A22 = -1.6009E-10
A23 = 7.994E-12  
A30 = 1.100E-10  
A31 = 6.651E-12  
A32 = -3.391E-13 
B00 = -1.922E-2  
B01 = -4.42E-5   
B10 = 7.3637E-5  
B11 = 1.7950E-7  
D00 = 1.727E-3   
D10 = -7.9836E-6


def c_chen(T,S,P):
    Cw = (C00 + C01*T + C02*T**2 + C03*T**3 + C04*T**4 + C05*T**5) + \
         (C10 + C11*T + C12*T**2 + C13*T**3 + C14*T**4)*P + \
  	 (C20 + C21*T + C22*T**2 + C23*T**3 + C24*T**4)*P**2 + \
  	 (C30 + C31*T + C32*T**2)*P**3
 
    A = (A00 + A01*T + A02*T**2 + A03*T**3 + A04*T**4) + \
  	(A10 + A11*T + A12*T**2 + A13*T**3 + A14*T**4)*P + \
  	(A20 + A21*T + A22*T**2 + A23*T**3)*P**2 + \
  	(A30 + A31*T + A32*T**2)*P**3
    
    B = B00 + B01*T + (B10 + B11*T)*P
 
    D = D00 + D10*P 

    return Cw + A*S + B * pow(S,3./2) + D*S**2

def c(T, S, Z, phi):
    c = 1402.5 + 5*T - 5.44e-2 * T**2 + 2.1e-4 * T**3 \
        + 1.33*S - 1.23e-2 * S*T + 8.7e-5*S*T**2 \
        + 1.56e-2*Z + 2.55e-7 *Z**2 - 7.3e-12*Z**3 \
        + 1.2e-6 *Z*(phi-45.) - 9.5e-13 * T*Z**3 \
        + 3e-7*T**2 *Z + 1.43e-5*S*Z
    
    return c


def main():
    
    new = False
    if new :
        C = c(13, 3.4, 2000., 45)
        a = np.arange(-0.05, 0.05, 0.01)
        b = np.arange(-10., 11., 1.)
        d = np.arange(-0.2, 0.21, 0.01)
        
        x = np.array( [13+i for i in a])
        y = np.array( [c(13+i, 3.4, 2000., 45) for i in a] )
        plt.subplot(311)
        plt.plot(x,y/C,'b')
        plt.xlabel("temperature [C]")
        plt.ylabel("c/c$_o$ ")
        plt.grid()
        
        x = np.array( [3.4+i for i in d] )
        y = np.array( [c(13, 3.4+i, 2000., 45) for i in d] )
        
        plt.subplot(312)
        plt.plot(x, y/C, 'r')
        plt.xlabel("salinity [%]")
        plt.ylabel("c/c$_o$ ")
        plt.grid()
        
        x = np.array( [2000+i for i in b] )
        y = np.array( [c(13, 3.4, 2000.+i, 45) for i in b] )
        
        plt.subplot(313)
        plt.plot(x, y/C, 'g')
        plt.xlabel("depth [m]")
        plt.ylabel("c/c$_o$ ")
        plt.grid()
        plt.tight_layout()
        
    else:
        C = c_chen(13, 34, 250.)
        print "speed of sound: ", C
        a = np.arange(-0.07, 0.07, 0.01)
        b = np.arange(-1., 1.1, .1)
        d = np.arange(-.2, .21, 0.001)
        
        x = np.array( [13+i for i in a])
        y = np.array( [c_chen(13+i, 34, 250.) for i in a] )
        plt.subplot(311)
        plt.plot(x,y/C,'b')
        plt.xlabel("temperature [C]")
        plt.ylabel("c/c$_o$ ")
        plt.grid()
        
        x = np.array( [34.+i for i in d] )
        y = np.array( [c_chen(13, 34+i, 250.) for i in d] )
        
        plt.subplot(312)
        plt.plot(x, y/C, 'r')
        plt.xlabel("salinity [psu]")
        plt.ylabel("c/c$_o$ ")
        plt.grid()

        x = np.array( [250+i for i in b] )
        y = np.array( [c_chen(13, 34, 250.+i) for i in b] )
        
        plt.subplot(313)
        plt.plot(x, y/C, 'g')
        plt.xlabel("pressure [bar]")
        plt.ylabel("c/c$_o$ ")
        plt.grid()
        plt.tight_layout()



#
#    x,y = np.mgrid[slice(2000, 2010, 1.),
#                    slice(3.3 , 3.5, 0.01)]
#    z = c(13, y, x, 34)
#    
#
#    fig = plt.figure()
#    ax = Axes3D(fig)
#    ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap=cm.coolwarm)
#    plt.ylabel("Depth [m]")
#    plt.xlabel("Salinity [%]")
#    ax.set_zlabel("Speed of sound")
    plt.show()
    
if __name__ == "__main__":
    main()

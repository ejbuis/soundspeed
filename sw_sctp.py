import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def sw_sal78(c, t, p):
    """
    Simplified version of the original SAL78 function from Fofonoff and Millard
    (1983). This does only the conversion from conductivity, temperature and
    pressure to salinity. Returns zero for conductivity values below 0.0005.
    Parameters
    ----------
    c : ndarray
        Conductivity (S m{-1})
    t : ndarray
        Temperature (degrees Celsius IPTS-68)
    p : ndarray
        Pressure (decibars)
    Returns
    -------
    s : salinity (PSU-78)
    Notes
    -----
    The Conversion from IPTS-68 to ITS90 is:
        T90 = 0.99976 * T68
        T68 = 1.00024 * T90
    These constants are defined here as c90 (0.99976) and c68 (1.00024).
    """

    p = p / 10

    C1535 = 1.0

    DT = t - 15.0

    # Convert conductivity to salinity
    rt35 = np.array((((1.0031E-9 * t - 6.9698E-7) * t + 1.104259E-4)
            * t + 2.00564E-2) * t + 0.6766097)
    a0 = np.array(-3.107E-3 * t + 0.4215)
    b0 = np.array((4.464E-4 * t + 3.426E-2) * t + 1.0)
    c0 = np.array(((3.989E-12 * p - 6.370E-8) * p + 2.070E-4) * p)

    R = np.array(c / C1535)
    RT = np.sqrt(np.abs(R / (rt35 * (1.0 + c0 / (b0 + a0 * R)))))

    s = np.array(((((2.7081 * RT - 7.0261) * RT + 14.0941) * RT + 25.3851)
        * RT - 0.1692) * RT + 0.0080 + (DT / (1.0 + 0.0162 * DT)) *
            (((((-0.0144 * RT + 0.0636) * RT - 0.0375) *
            RT - 0.0066) * RT - 0.0056) * RT + 0.0005))

    # Zero salinity trap
    if len(s.shape) > 0:
        s[c < 5e-4] = 0

    return s


def main():
#    Salinity=40.00000 for CND=1.888091, T=40 DEG C, P=10000 DECIBAR
#    sw_sal78(1.888091, 40, 10000)

    c = 1.0
    t = 13
    p = 2500
    S = sw_sal78(c, t, p)
    print("Salinity ", S)
    
#    a = np.arange(-0.005, 0.005, 0.0001)
    a = np.arange(-0.15, 0.5, 0.1)
    b = np.arange(-0.3, 0.3, 0.001)
    d = np.arange(-10., 10., 1.)
    
    x = np.array( [c+i for i in a])
    y = np.array( [sw_sal78(c+i, t, p) for i in a] )

    X = [x[0], x[-1]]
    Yu = [1.003, 1.003]
    Yl = [0.997, .997]    
    ax = plt.subplot(311)
    ax.add_line( Line2D(X, Yu, lw = 1.5, color='r', linestyle = '--') )
    ax.add_line( Line2D(X, Yl, lw = 1.5, color='r', linestyle = '--') )
    plt.plot(x,y/S, 'b')
    plt.xlabel("conductivity [S/m]")
    plt.ylabel("S/S$_o$ ")
    plt.grid()
    
    
    x = np.array( [t+i for i in b] )
    y = np.array( [sw_sal78(c, t+i, p) for i in b] )
    X = [x[0], x[-1]]
    
    ax = plt.subplot(312)
#    ax.add_line( Line2D(X, Yu, lw = 2.5, color='r', linestyle = ':') )
#    ax.add_line( Line2D(X, Yl, lw = 2.5, color='r', linestyle = ':') )
    plt.plot(x, y/S, 'r')
    plt.xlabel("temperature [$^o$C]")
    plt.ylabel("S/S$_o$ ")
    plt.grid()
    
    x = np.array( [p+i for i in d] )
    y = np.array( [sw_sal78(c, t, p+i) for i in d] )
    X = [x[0], x[-1]]
    
    ax = plt.subplot(313)
#    ax.add_line( Line2D(X, Yu, lw = 2.5, color='r', linestyle = ':') )
#    ax.add_line( Line2D(X, Yl, lw = 2.5, color='r', linestyle = ':') )
    plt.plot(x/10, y/S, 'g')
    plt.xlabel("pressure [bar]")
    plt.ylabel("S/S$_o$ ")
    plt.grid()
    plt.tight_layout()
        
    plt.show()



if __name__ == "__main__":
    main()

    

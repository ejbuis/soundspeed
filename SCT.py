import numpy as np

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
    c = 1.888091
    t = 40
    p = 10000
    print sw_sal78(c, t, p)



if __name__ == "__main__":
    main()

    

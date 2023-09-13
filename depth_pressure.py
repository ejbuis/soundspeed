import math
import numpy as np
import sys

# Function Cond2Sal converts a conductivity value of seawater to a value
# of the pratical-salinity-scale 1978 (PSS-78) for given values of
# conductivity, temperature and pressure. Result is returned as
# parameter in aSalinity. A returned boolean result TRUE of the
# function indicates that the result is reliable.
# UNITS:
#   PRESSURE      Press          DECIBARS
#   TEMPERATURE   Temp           DEG CELSIUS IPTS-68
#   CONDUCTIVITY  aConductivity  S/m
#   SALINITY      aSalinity      PSS-78
# ----------------------------------------------------------
# CHECKVALUES:
#   2.) aSalinity=40.00000 for CND=1.888091, T=40 DEG C, P=10000 DECIBARS
# ----------------------------------------------------------
# SAL78 RATIO: RETURNS ZERO FOR CONDUCTIVITY RATIO: < 0.0005
# ----------------------------------------------------------
# This source code is based on the original fortran code in:
#   UNESCO technical papers in marine science 44 (1983) -
#   'Algorithms for computation of fundamental properties of seawater'
# ----------------------------------------------------------
# Written in object pascal by:
#   Dr. Jan Schulz, 26. May 2008, www.code10.info
#  convert to python, Ernst-Jan Buis

def SAL(XR, XT):
    # PRACTICAL SALINITY SCALE 1978 DEFINITION WITH TEMPERATURE
    # CORRECTION;XT :=T-15.0; XR:=SQRT(RT);
    
    return  ((((2.7081*XR-7.0261)*XR+14.0941)*XR+25.3851)*XR \
             - 0.1692)*XR+0.0080 + (XT/(1.0+0.0162*XT))* \
             (((((-0.0144*XR + \
                  0.0636)*XR-0.0375)*XR-0.0066)*XR-0.0056)*XR+0.0005)
    

def RT35(XT): 
    # FUNCTION RT35: C(35,T,0)/C(35,15,0) VARIATION WITH TEMPERATURE
    return (((1.0031E-9 * XT - 6.9698E-7) * XT + 1.104259E-4) * XT
            + 2.00564E-2) * XT + 0.6766097

def C(XP):
    # C(XP) POLYNOMIAL CORRESPONDS TO A1-A3 CONSTANTS: LEWIS 1980
    return ((3.989E-15 * XP - 6.370E-10) * XP + 2.070E-5) * XP

def B(XT): 
    return (4.464E-4 * XT + 3.426E-2) * XT + 1.0;
    

def A(XT):
    # A(XT) POLYNOMIAL CORRESPONDS TO B3 AND B4 CONSTANTS: LEWIS 1980
    return -3.107E-3 * XT + 0.4215;

def Cond2Sal78(aConductivity, Temp, Press):
    # start conversion
    DT = Temp - 15
    aSalinity = aConductivity/4.2914
    RT = aSalinity / (RT35 (Temp) * (1.0 + C(Press) /
                                     (B(Temp) + A(Temp) * aSalinity)))
    RT = math.sqrt(np.abs(RT))
    aSalinity = SAL(RT, DT)

    # control, whether result is in the validity range of PSS-78
    if (aSalinity < 2) or (aSalinity > 42):
        sys.exit("wrong values!! No result")
    return aSalinity

def main():
    #
    aConductivity = 1.888091
    Temp = 25
    Press = 10000
    print Cond2Sal78(aConductivity, Temp, Press)
    
    
if __name__ == "__main__":
    main()


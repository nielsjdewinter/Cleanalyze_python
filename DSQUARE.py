import numpy as np

def DSQUARE(Yin, WINDOW):
    """
    Yin           The data that is to be analysed
    WINDOW        Size of the window for which mean and variance will be
                  calculated before and after the point of interest
    """
    DSq = np.full_like(Yin, np.nan, dtype=np.float64)
    for i in range(WINDOW, len(Yin) - WINDOW):
        MUb = np.mean(Yin[i-WINDOW:i])  # MUb is average from window before point i
        VARb = np.var(Yin[i-WINDOW:i])
        MUa = np.mean(Yin[i+1:i+1+WINDOW])  # MUa is average from window after point i
        VARa = np.var(Yin[i+1:i+1+WINDOW])
        DSq[i] = (MUb - MUa)**2 / (VARb + VARa)
    return DSq

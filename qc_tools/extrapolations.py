import numpy as np
import scipy.stats as st

def X3(X, Y):
    if len(X) != len(Y):
        raise Exception('The length of X must be equal to that of Y')

    Xm3 = np.array([float(x)**(-3) for x in X])
    Yp = np.array(Y)
    slope, intercept, r, p, s = st.linregress(Xm3, Yp)

    return(intercept)

def exp(X, Y, alpha = 1.63):
    if len(X) != len(Y):
        raise Exception('The length of X must be equal to that of Y')

    Xexp = np.exp(-alpha * np.array(X))
    Yp = np.array(Y)

    slope, intercept, r, p, s = st.linregress(Xexp, Yp)

    return(intercept)

if __name__ == '__main__':
    X = [3, 4]
    Y = [-31.7977581,     -31.850643]

    print(X3(X, Y))
    print(exp(X, Y))

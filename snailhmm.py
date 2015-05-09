import scipy
import numpy as np

numIntervals = 20

tmax = 10.0

t = []

for i in range(numIntervals):
    t.append(0.1*(np.exp(float(i)/numIntervals*np.log(1+10*tmax))-1.0))
t.append(tmax)

t = np.array(t)

def pi_k(k, t):
    t1 = t[k]
    t2 = t[k+1]
    pik = 1.0/8.0*(np.exp(-3.0*t2)-np.exp(-3.0*t1)+np.exp(-t1)*(9.0+6.0*t1)-3.0*np.exp(-t2)*(3.0+2.0*t2))
    return pik

def sigma_k(k, t):
    t1 = t[k]
    t2 = t[k+1]
    sigmak = np.exp(-t1)-np.exp(-t2)
    return sigmak

def qkl(k,l,t):
    tk1 = t[k]
    tk2 = t[k+1]
    tl1 = t[l]
    tl2 = t[l+1]
    if k < l:
        qkl_integral = 3.0/8.0*np.exp(-2.0*tk1-2.0*tk2-tl1-tl2)*(np.exp(tl1)-np.exp(tl2))*(np.exp(2*tk2)-np.exp(2*tk1)+2.0*np.exp(2.0*(tk1+tk2))*(tk1-tk2))
        qkl = qkl_integral/pi_k(k,t)
        return qkl
    if k > l:
        qkl_integral = 3.0/8.0*np.exp(-2.0*tl1-2.0*tl2-tk1-tk2)*(np.exp(tk1)-np.exp(tk2))*(np.exp(2*tl2)-np.exp(2*tl1)+2.0*np.exp(2.0*(tl1+tl2))*(tl1-tl2))
        qkl = qkl_integral/pi_k(k,t)
        return qkl
    if k == l:
        t1 = tk1
        t2 = tk2
        qkl_integral = 1./4.*np.exp(-3.0*(t1+t2))*(-np.exp(3*t1)+2.0*np.exp(3.0*t2)*(3.0*np.exp(2.0*t1)-1)+3.0*np.exp(t1+2.0*t2)*(1.0+2.0*np.exp(2*t1)*(-1+t1-t2)))
        qkl = qkl_integral/pi_k(k,t)
        return qkl

def qkl2(k,l,t):
    tk1 = t[k]
    tk2 = t[k+1]
    tl1 = t[l]
    tl2 = t[l+1]
    if k < l:
        qkl_integral = 3.0/8.0*np.exp(-2.0*tk1-2.0*tk2-tl1-tl2)*(np.exp(tl1)-np.exp(tl2))*(np.exp(2*tk2)-np.exp(2*tk1)+2.0*np.exp(2.0*(tk1+tk2))*(tk1-tk2))
        qkl = qkl_integral/pi_k(k,t)
        return qkl
    if k > l:
        qkl_integral = 3.0/8.0*np.exp(-2.0*tl1-2.0*tl2-tk1-tk2)*(np.exp(tk1)-np.exp(tk2))*(np.exp(2*tl2)-np.exp(2*tl1)+2.0*np.exp(2.0*(tl1+tl2))*(tl1-tl2))
        qkl = qkl_integral/pi_k(k,t)
        return qkl
    if k == l:
        t1 = tk1
        t2 = tk2
        qkl_integral = 1./4.*np.exp(-3.0*(t1+t2))*(-np.exp(3.0*t1)+2.0*np.exp(3.0*t2)*(3.0*np.exp(2.0*t1)-1.0)+3*np.exp(t1+2.0*t2)*(1.0+2.0*np.exp(2.0*t1)*(t1-t2-1)))
        qkl = qkl_integral/pi_k(k,t)
        return qkl

[(sum([qkl(j,l,t) for l in range(numIntervals)]),sum([qkl(j,l,t) for l in range(numIntervals)])) for j in range(20)]
# fix/check this...
sum([qkl(19,l,t) for l in range(numIntervals)])
sum([qkl2(19,l,t) for l in range(numIntervals)])
[qkl(19,l,t) for l in range(numIntervals)]


[pi_k(k, t) for k in range(numIntervals)]

[sigma_k(k, t) for k in range(numIntervals)]



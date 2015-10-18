import scipy
import scipy.integrate
import scipy.optimize
import numpy as np
import sys

numIntervals = 50

tmax = 15.0

t = []

for i in range(numIntervals-1):
    t.append(0.1*(np.exp(float(i)/(numIntervals-1)*np.log(1+10*tmax))-1.0))
t.append(tmax)
t.append(float('inf'))

t = np.array(t)

def pi_k(k, t):
    if k >= len(t)-1:
        return None
    t1 = t[k]
    t2 = t[k+1]
    if t2 != float('inf'):
        pik = 1.0/8.0*(np.exp(-3.0*t2)-np.exp(-3.0*t1)+np.exp(-t1)*(9.0+6.0*t1)-3.0*np.exp(-t2)*(3.0+2.0*t2))
        return pik
    elif t2 == float('inf'):
        pik = 1.0/8.0*np.exp(-3.0*t1)*(-1.0+np.exp(2.0*t1)*(9.0+6.0*t1))
        return pik

def sigma_k(k, t):
    if k >= len(t)-1:
        return None
    t1 = t[k]
    t2 = t[k+1]
    sigmak = np.exp(-t1)-np.exp(-t2)
    return sigmak

def qkl(k,l,t):
    if k >= len(t)-1 or l >= len(t)-1:
        return None
    tk1 = t[k]
    tl1 = t[l]
    tk2 = t[k+1]
    tl2 = t[l+1]
    if k < l:
        if tl2 != float('inf'):
            qkl_integral = 3.0/8.0*np.exp(-2.0*tk1-2.0*tk2-tl1-tl2)*(np.exp(tl1)-np.exp(tl2))*(np.exp(2*tk2)-np.exp(2*tk1)+2.0*np.exp(2.0*(tk1+tk2))*(tk1-tk2))
            qkl = qkl_integral/pi_k(k,t)
            return qkl
        elif tl2 == float('inf'):
            qkl_integral = 3.0/8.0 * np.exp(-tl1)*(-np.exp(-2.0*tk1)+np.exp(-2.0*tk2)-2.0*tk1+2.0*tk2)
            qkl = qkl_integral/pi_k(k,t)
            return qkl
    if k > l:
        if tk2 < float('inf'):
            qkl_integral = 3.0/8.0*np.exp(-2.0*tl1-2.0*tl2-tk1-tk2)*(np.exp(tk1)-np.exp(tk2))*(np.exp(2*tl2)-np.exp(2*tl1)+2.0*np.exp(2.0*(tl1+tl2))*(tl1-tl2))
            qkl = qkl_integral/pi_k(k,t)
            return qkl
        elif tk2 == float('inf'):
            qkl_integral = 3.0/8*np.exp(-tk1)*(-np.exp(-2.0*tl1)+np.exp(-2.0*tl2)-2.0*tl1+2*tl2)
            qkl = qkl_integral/pi_k(k,t)
            return qkl
    if k == l:
        t1 = tk1
        t2 = tk2
        if t2 != float('inf'):
            qkl_integral = 1./4.*np.exp(-3.0*(t1+t2))*(-np.exp(3*t1)+2.0*np.exp(3.0*t2)*(3.0*np.exp(2.0*t1)-1)+3.0*np.exp(t1+2.0*t2)*(1.0+2.0*np.exp(2*t1)*(-1+t1-t2)))
            qkl = qkl_integral/pi_k(k,t)
            return qkl
        elif t2 == float('inf'):
            qkl_integral = 1./2.*np.exp(-3.*t1)*(-1.+3.*np.exp(2.*t1))
            qkl = qkl_integral/pi_k(k,t)
            return qkl

#print [sum([qkl(j,l,t) for l in range(numIntervals)]) for j in range(numIntervals)]
# beautiful, working again

def lambda_(s):
    return 1./4.*(1.-np.exp(-2.*s)+2.*s)

def pkLTl_integrand(t,s,rho):
    return (1.0-np.exp(-lambda_(s)*rho))*(2.0*np.exp(s-t)*(1.0-np.exp(-2.0*s)))/(1.0-np.exp(-2.*s)+2.*s)*np.exp(-s)

def pkLTl(k,l,t,rho):
    if k >= len(t)-1 or l >= len(t)-1:
        return None
    tk1 = t[k]
    tl1 = t[l]
    tk2 = t[k+1]
    tl2 = t[l+1]
    integral, err = scipy.integrate.dblquad(pkLTl_integrand, tk1,tk2, lambda x: tl1, lambda x: tl2, args=(rho,))
    pkl = integral/sigma_k(k,t)
    return pkl

def pkGTl_integrand(t,s,rho):
    return (1.0-np.exp(-lambda_(s)*rho))*(2.0*(1.0-np.exp(-2.0*t)))/(1.0-np.exp(-2.*s)+2.*s)*np.exp(-s)

def pkGTl(k,l,t,rho):
    if k >= len(t)-1 or l >= len(t)-1:
        return None
    tk1 = t[k]
    tl1 = t[l]
    tk2 = t[k+1]
    tl2 = t[l+1]
    integral, err = scipy.integrate.dblquad(pkGTl_integrand, tk1,tk2, lambda x: tl1, lambda x: tl2, args=(rho,))
    pkl = integral/sigma_k(k,t)
    return pkl

def pkEQl_integrand(s,rho):
    return np.exp(-lambda_(s)*rho)*np.exp(-s)

def pkEQl(k, l, t, rho):
    if k != l:
        return None
    if k >= len(t)-1 or l >= len(t)-1:
        return None
    t1 = t[k]
    t2 = t[k+1]

    integral = 0.0

    integral0, err0 = scipy.integrate.dblquad(pkGTl_integrand, t1,t2, lambda x: t1, lambda x: x, args=(rho,))

    integral += integral0

    integral1, err1 = scipy.integrate.dblquad(pkLTl_integrand, t1,t2, lambda x: x, lambda x: t2, args=(rho,))
    
    integral += integral1

    integral2, err2 = scipy.integrate.quad(pkEQl_integrand, t1,t2, args=(rho,))
    
    integral += integral2

    return integral/sigma_k(k,t)


#def pkEQl(k, l, t, rho):

def pkl(k, l, t, rho):
    if k < l:
        return pkLTl(k,l,t,rho)
    if k > l:
        return pkGTl(k,l,t,rho)
    if k == l:
        return pkEQl(k,l,t,rho)

#print [sum([pkl(j,l,t, 0.25) for l in range(numIntervals)]) for j in range(numIntervals)]
# score

def ek1(k, x, t, theta, P):
    ek1 = 0
    if k >= len(t)-1:
        return None
    tk1 = t[k]
    tk2 = t[k+1]
    
    ek1 = np.exp(-tk1)-np.exp(-tk2)

    ek1 -= np.exp(-P*theta-tk1*(1+theta))/(1.0+theta)

    ek1 += np.exp(-P*theta-tk2*(1+theta))/(1.0+theta)

    if x == 1:
        return ek1
    elif x == 0:
        return 1-ek1

examplePositions = []

examplePosIn = open('example_data', 'r')
for line in examplePosIn:
    line = line.strip()
    examplePositions.append(float(line))

examplePositions = np.array(examplePositions)

# Li and Durbin uses rho ~ 0.01 per bin
# (1 cM per Mb, rho = 0.01*10000 per 1000000, which is 0.01 per 100 bp)

# example data has rho = 1000, so 100000 bins
# has theta = 1000
# and P = 0.2

intervals = np.linspace(0.0, 1.0, num=100000)

data = []

for start, end in zip(intervals[:-1], intervals[1:]):
    found = False
    for pos in examplePositions:
        if start < pos and pos < end:
            found = True
            break
    if found:
        data.append(1)
    else:
        data.append(0)

data = np.array(data)

###########################
# inference

# twice as high as they should be...
originalTheta = 0.2
originalRho = 0.2
originalP = 0.5

#pkl_mat = np.matrix([[pkl(k,l,t,originalRho) for l in range(numIntervals)] for k in range(numIntervals)])

# calculate first alphas

def get_likelihood(args, x, t):

    x = data

    theta = args[0]
    rho = args[1]
    P = args[2]

    numSites = len(x)

    pkl_mat = [[pkl(i,j,t,rho) for j in range(numIntervals)] for i in range(numIntervals)]
    for row in pkl_mat:
        rowSum = 0.0
        for entry in row:
            rowSum += entry
        for j in range(len(row)):
            row[j] /= rowSum

    logsums = []

    alpha0 = np.array([sigma_k(k,t)*ek1(k, x[0], t,theta,P) for k in range(numIntervals)])
    alphasum = sum(alpha0)
    logsums.append(-np.log(alphasum))
    alpha0 /= alphasum
    curAlphas = alpha0
    for k in range(1, numSites):
        alphaPrime = np.array([sum([curAlphas[i]*pkl_mat[i][j] for i in range(numIntervals)])*ek1(j,x[k], t, theta, P) for j in range(numIntervals)])
        alphasum = sum(alphaPrime)
        logsums.append(-np.log(alphasum))
        alphaPrime /= alphasum
        curAlphas = alphaPrime

    loglike = -1.0*sum(logsums)

    return loglike

def neg_like(args, x ,t):
    return -1.0*get_likelihood(args, x, t)

#max_like = scipy.optimize.fmin_powell(neg_like, (originalTheta, originalRho, originalP), args = (data, t), full_output = True, maxfun = 20)

bounds = [(1e-12, None)] * 3
max_like = scipy.optimize.fmin_tnc(neg_like, (originalTheta, originalRho, originalP), args = (data, t), bounds = bounds, approx_grad = True, disp=5)
#max_like = scipy.optimize.fmin_l_bfgs_b(neg_like, (originalTheta, originalRho, originalP), args = (data, t), approx_grad = True, bounds = bounds, iprint = 0)

print max_like



# Homework 3 Numerical Method
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#                                                             Number 1
# Solve the following integral using numerical methods presented here:
# Integrate(1+x**2), first we define function


def fx(x):  # call function 'remember fx FOR SCALARS, fxx for vectors'
    y = math.sqrt(1+x**2)  # define function
    return y


# define limits
a = 0  # lower bound
b = 1  # upper bound

seg = 100  # break into segments
dx = (b-a)/seg  # x length of each segment
x = np.arange(a+dx, b, dx)  # make all the x segments not including a and b
fxx = np.vectorize(fx)  # since we're passing vector we cant use fx which...
# ...expects a scalar so we use a a vectorized function version instead

#                                   1a. using ➢ Trapezoidal Rule(Composite)

sum_trapz = np.sum(fxx(x))  # sum of the y(vlaues of x internal segments)
trapz = dx/2*(fx(a)+2*sum_trapz+fx(b))  # compute integral using trapezoid
print(trapz)  # 1.148 is answer using 1000 broken segment

#                                               1b. using ➢ Simpson’s rule

twos_simsum = 2*(np.sum(fxx(x[2:-1:2])))  # sum for all the 'EVEN' internal y segments
fours_simpsum = 4*(np.sum(fxx(x[1:-1:2])))  # sum for all the 'ODD' internal y segments
simp = (dx/3)*(fx(a)+twos_simsum+fours_simpsum+fx(b))  # compute integral using simpsomn.
print(simp)  # 1.145 is answer using 1000 broken segment

#                                               1c. using ➢ Romberg Integration
unit = 10  # how many rows for romberg table
romb = []  # the list must be defined before elements can be added
h = (b-a)  # caluculate the dx
romb.append([h/2*(fxx(a)+fxx(b))])  # start with first index in romb
for i in range(1, unit+1):
    h = h/2
    sum = 0
    for k in range(1, 2**i, 2):
        sum += fxx(a+k*h)

    rowi = [0.5*romb[i-1][0] + sum*h]  # start the next row fr iteration
    for j in range(1, i+1):  # compute iteration for next row
        rij = rowi[j-1] + (rowi[j-1]-romb[i-1][j-1])/(4**j-1)
        rowi.append(rij)  # Append to rowi
    romb.append(rowi)  # finally append to romb

print(romb[9][1])
#
#
#                                          1d. using ➢ Gauss Quadrature(2 point)


def gx(x):  # define function
    z = np.sqrt(x**2+1)
    return z


gauss = []
gxx = np.vectorize(gx)

# define weights coordinates from gauss quadrature 2 point formula
w_1 = w_2 = 1
z_1 = 1/(3**0.5)
z_2 = -z_1
x_1 = (b-a)/2*z_1+(b+a)/2
x_2 = (b-a)/2*z_2+(b+a)/2
gauss = (b-a)/2*(w_1*gxx(x_1)+w_2*gxx(x_2))
print(gauss)  # 1.1478

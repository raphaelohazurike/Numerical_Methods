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
#----------------------------------------------------------------------------------------------------
# Homework 3 Numerical Method
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# NUMBER 2.
# pump force is directly proportional to water distance travelled
#       Compute the work done by the pump moving water over a distance of 10 m.
# F = 2 x1.5 (where x is in meters and F is in kN)#
#                                                   2a. Compare Gauss Quadrature


def gp(x):  # define function
    y = 2*x**1.5
    return y


gauss = []  # initialize answer
gpx = np.vectorize(gp)
a = 0  # lower limits
b = 10  # upper limit
# define weights coordinates from gauss quadrature 2 point formula
w_1 = w_2 = 1
z_1 = 1/(3**0.5)
z_2 = -z_1
x_1 = (b-a)/2*z_1+(b+a)/2
x_2 = (b-a)/2*z_2+(b+a)/2
gauss = (b-a)/2*(w_1*gpx(x_1)+w_2*gpx(x_2))  # quad formula
print(gauss)  # 252.20kN
#                                                   # 2b. Trapezoidal Rule.
# Integrate(2x^1.5), first we define function


def px(x):  # call function 'remember fx FOR SCALARS, fxx for vectors'
    y = 2*x**1.5  # define function
    return y


# define limits
a = 0  # lower bound
b = 10  # upper bound

seg = 100  # break into segments
dx = (b-a)/seg  # x length of each segment
x = np.arange(a+dx, b, dx)  # make all the x segments not including a and b
pxx = np.vectorize(px)  # since we're passing vector we cant use fx which...
# ...expects a scalar so we use a a vectorized function version instead

sum_trapz = np.sum(pxx(x))  # sum of the y(vlaues of x internal segments)
wd = dx/2*(px(a)+2*sum_trapz+px(b))  # compute integral using trapezoid
print(wd)  # 252.9kN is answer using 100 broken segment
#                           2c. COMPARISON BETWEEN tRAPEZOID AND GAUSS QUADRATURE
# Gauss Quadrature has just 2 point, and is almost accurate, another iteration of
# of 10 or more poitns would give a more correct answer than the trapezoidal rule
#--------------------------------------------------------------------------------------------------------------------
# Homework 3 Numerical Method
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# NUMBER 3.
# Plate is 8 m long and 6 m wide.
# Compute average temperature


def Txy(x, y):  # define function
    z = 2*x*y + 2*x - x**2 - 2*y**2 + 72
    return z


#
tv = np.vectorize(Txy)  # vectorize funtion
# # # integral limits
xa = 0
xb = 8
ya = 0
yb = 6
seg = 1000
dx = (xb-xa)/seg
dy = (yb-ya)/seg
X = np.arange(xa+dx, xb, dx)  # make x array
Y = np.arange(ya+dx, yb, dx)  # make y array
X, Y = np.meshgrid(X, Y)
XX = X.flatten()  # make 1D
YY = Y.flatten()

# integrate
z = tv(XX, YY)  # calculate temperature values
tz = (dx*dy)/2*(np.sum(z))  # trapezoidal rule
print(tz)  # 1055
#
#                                    3b Solve using Romberg Integration(25 Points)
unit = 10  # how many rows for romberg table

romb = []  # the list must be defined before elements can be added

romb.append([dx/2*(tv(XX, YY))])  # start with first index in romb
for i in range(1, unit+1):
    dx = dx/2
    sum = 0
    for k in range(1, 2**i, 2):
        sum += fxx(a+k*dx)

    rowi = [0.5*romb[i-1][0] + sum*h]  # start the next row fr iteration
    for j in range(1, i+1):  # compute iteration for next row
        rij = rowi[j-1] + (rowi[j-1]-romb[i-1][j-1])/(4**j-1)
        rowi.append(rij)  # Append to rowi
    romb.append(rowi)  # finally append to romb

print(romb[9][1])
#                       3c Solve using Trapezoidal Rule + Gaussian Quadrature(25 points)
gauss = []  # initialize answer
# define weights coordinates from gauss quadrature 2 point formula
YY = np.arange(ya, yb+dx, dx)  # create Y array
Ny = len(YY)
np.zeros(Ny)
w_1 = w_2 = 1
z_1 = 1/(3**0.5)
z_2 = -z_1
x_1 = (b-a)/2*z_1+(b+a)/2
x_2 = (b-a)/2*z_2+(b+a)/2
gauss = (b-a)/2*(w_1*gpx(x_1)+w_2*gpx(x_2))  # quad formula
print(gauss)  #


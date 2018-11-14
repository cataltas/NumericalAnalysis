# Project on polynomial interpolation
# Written by Camille Taltas
from numpy import *
from  matplotlib.pyplot import *

#Defining functions for exercises 1 and 2
def g(x):
    return sin(2.0*pi*x)/x

def f(x):
    return (1-x**3)*sin(2.0*pi*x)

def g_limit(x):
    return 2.0*pi*cos(2.0*pi*x)


def polynomial_interpolation(func, interval, degree,limit_func=0):

    #Initializing Step Size
    step_size = (interval[1]-interval[0])/(degree)
    data_points = []
    data_points_f = []
    i = interval[0]
    #Finding all the data points and corresponding value and applying l'Hopital for undefined values 
    while i<= interval[1]:
        data_points.append(i)
        if isnan(func(i)) == True:
            data_points_f.append(limit_func(i))
        else:
            data_points_f.append(func(i))
        i = i+step_size
    #Finding the corresponding Vandermonde matrix to solve the system of equations and find the a's. 
    A = vander(data_points,increasing=True)
    B = data_points_f
    a_values = linalg.solve(A,B)
    #Finding the data points of our polynomial to then plot the graph of both interpolations 
    xvals = []
    yvals = []
    i = interval[0]
    while i<=interval[1]:
        y_val = 0
        for j in xrange (0,len(a_values)):
            y_val = y_val + a_values[j]*(i**j)
        yvals.append(y_val)
        xvals.append(i)
        i = i + 0.001
    plot(xvals,yvals, label=degree)
    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    title("Polynomial Interpolation")
    return yvals

#Question 1 Part 1 Graphing 
polynomial_interpolation(g,[-1,1.5],5,g_limit)
polynomial_interpolation(g,[-1,1.5],10,g_limit)
show()

#Question 1 Part 2 Graphing
polynomial_interpolation(f,[-1.0,2.0],6.0)
polynomial_interpolation(f,[-1.0,2.0],10.0)
#Plot f(x)
xvals = []
yvals = []
i = -1
while i<=2:
    y_val = f(i)
    yvals.append(y_val)
    xvals.append(i)
    i = i + 0.001
plot(xvals,yvals, label = "f(x)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#Show interpolants and graph
show()

def polynomial_interpolation_chebyshev(func,interval,degree):
    #Initializing the data points based on Chebyshev roots 
    data_points = []
    data_points_f =[]
    for i in xrange(1,(degree+2)):
        x = (interval[0]+interval[1])/2.0 + (interval[1]-interval[0])/2.0*cos(((2.0*i-1.0)/(2.0*(degree+1))*pi))
        f_x = func(x)
        data_points.append(x)
        data_points_f.append(f_x)
    #Finding the corresponding Vandermonde matrix to solve the system of equations and find the a's.
    A = vander(data_points,increasing=True)
    B = data_points_f
    a_values = linalg.solve(A,B)
    #Finding the data points of our polynomial to then plot the graph of both interpolations 
    xvals = []
    yvals = []
    i = interval[0]
    while i<=interval[1]:
        y_val = 0
        for j in xrange (0,len(a_values)):
            y_val = y_val + a_values[j]*(i**j)
        yvals.append(y_val)
        xvals.append(i)
        i = i + 0.001
    plot(xvals,yvals, label = degree)
    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    title("Polynomial Interpolation with Chebyshev nodes")
    return yvals

#Question 2 Part 1 
polynomial_interpolation_chebyshev(g,[-1,1.5],5)
polynomial_interpolation_chebyshev(g,[-1,1.5],10)
show()
#Question 2 Part 2
polynomial_interpolation_chebyshev(f,[-1.0,2.0],6)
polynomial_interpolation_chebyshev(f,[-1.0,2.0],10)
#Plot f(x)
xvals = []
yvals = []
i = -1
while i<=2:
    y_val = f(i)
    yvals.append(y_val)
    xvals.append(i)
    i = i + 0.001
plot(xvals,yvals, label = "f(x)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#Show interpolants and graph
show()

#Defining b_0 to calculate S_3
def b_0(x,a,h):
    if x-a <= -2*h:
        return 0
    elif x-a<= -h:
        return 1.0/6.0*((2*h+(x-a))**3)
    elif x-a<=0:
        return 2.0/3.0*(h**3)-1.0/2.0*((x-a)**2)*(2*h+(x-a))
    elif x-a<=h:
        return 2.0/3.0*(h**3)-1.0/2.0*((x-a)**2)*(2*h-(x-a))
    elif x-a<=2*h:
        return 1.0/6.0*((2*h-(x-a))**3)
    elif x-a>=2.0*h:
        return 0

def natural_cubic_spline(func, interval,n_points,limit_func=0):
    #Initializing Step Size 
    step_size = (interval[1]-interval[0])/(n_points)
    data_points = []
    data_points_f = []
    i = interval[0]
    #Finding all the (n+1) data points and corresponding value and applying l'Hopital for undefined values 
    while i<= interval[1]:
        data_points.append(i)
        if isnan(func(i)) == True:
            data_points_f.append(limit_func(i))
        else:
            data_points_f.append(func(i))
        i = i+step_size
    #Setting up the matrix needed to solve the system of equations
    A = zeros((n_points+1,n_points+1))
    for i in xrange(0,n_points+1):
        for j in xrange (0,n_points+1):
            if j == 0:
                A[0,0]=1
            elif j==n_points:
                A[n_points,n_points]=1
            else:
                A[j,j-1]=1
                A[j,j]=4
                A[j,j+1]=1

    B=[]
    for j in xrange (0,n_points+1):
        if j == 0 | j == n_points:
            B.append(1/(step_size**3)*data_points_f[j])
        else:
            B.append(6/(step_size**3)*data_points_f[j])
    #Solving the system for the a_0 to a_n values and setting up the a_1, a_n+1 values assuming we have a natural spline
    a_values = linalg.solve(A,B)
    a_1 = 2*a_values[0]-a_values[1]
    a_n1 = 2*a_values[n_points]-a_values[n_points-1]
    #Finding the data points of our polynomial to then plot the graph of both interpolations using b_0
    xvals = []
    yvals = []
    i = interval[0]
    while i<=interval[1]:
        yval=0
        for k in xrange(-1,n_points+2):
            if k == -1:
                yval = yval + a_1*b_0(i-k*step_size,interval[0],step_size)
            elif k == n_points+1:
                yval = yval + a_n1*b_0(i-k*step_size, interval[0],step_size)
            else:
                yval = yval + a_values[k]*b_0(i-k*step_size,interval[0],step_size)
        yvals.append(yval)
        xvals.append(i)
        i = i + 0.001
    plot(xvals,yvals, label = n_points)
    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    title("Natural Spline")
    return yvals

#Question 3 Part 1
natural_cubic_spline(g,[-1.0,1.5],5,g_limit)
natural_cubic_spline(g,[-1.0,1.5],10,g_limit)
show()

#Question 3 Part 2
natural_cubic_spline(f,[-1.0,2.0],6)
natural_cubic_spline(f,[-1.0,2.0],10)
#Plot f(x)
xvals = []
yvals = []
i = -1
while i<=2:
    y_val = f(i)
    yvals.append(y_val)
    xvals.append(i)
    i = i + 0.001
plot(xvals,yvals, label = "f(x)")
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#Show interpolants and graph
show()

#Computing the actual maximum error 
def max_error(interpolant,func, interval):
    i = interval[0]
    yvals = []
    while i<=interval[1]:
        yvals.append(f(i))
        i = i + 0.001
    error = abs(array(yvals)-array(interpolant))
    return max(error)

#Computing the theoretical max error for the interpolating polyomial
#The maximum of the n+1th derivative is computed on Wolfram Alpha
def theoretical_error_inter(interval,degree, max_derivative):
    i = interval[0]
    max_product = 0
    #Finding the maximum product
    while i<=interval[1]:
        product = 1.0
        step_size = (interval[1]-interval[0])/(degree)
        k = interval[0]
        while k<= interval[1]:
            product = product * abs(i-k)
            k = k+step_size
        if product > max_product:
            max_product=product 
        i = i+0.001
    return max_product/math.factorial(degree+1)*max_derivative

#Computing the theoretical maximum error for the interpolating polynomial with Chebyshev nodes
#The maximum of the n+1th derivative is computed on Wolfram Alpha
def theoretical_error_chebyshev(degree, max_derivative):
    return 1.0/((2.0**degree)*math.factorial(degree+1))*max_derivative

print "The actual maximum error for the interpolating polynomial of degree 6 is", max_error(polynomial_interpolation(f,[-1.0,2.0],6),f,[-1.0,2.0])
print "The theoretical maximum error for the interpolating polynomial of degree 6 is", theoretical_error_inter([-1.0,2.0],6,3.78707*(10**6))

print "The actual maximum error for the interpolating polynomial of degree 10 is", max_error(polynomial_interpolation(f,[-1.0,2.0],10),f,[-1.0,2.0])
print "The theoretical maximum error for the interpolating polynomial of degree 10 is", theoretical_error_inter([-1.0,2.0],10,1.04791*(10**10))

print "The actual maximum error for the Chebyshev interpolating polynomial of degree 6 is", max_error(polynomial_interpolation_chebyshev(f,[-1.0,2.0],6),f,[-1.0,2.0])
print "The theoretical maximum error for the Chebyshev interpolating polynomial of degree 6 is", theoretical_error_chebyshev(6,3.78707*(10**6))

print "The actual maximum error for the Chebyshev interpolating polynomial of degree 10 is", max_error(polynomial_interpolation_chebyshev(f,[-1.0,2.0],10),f,[-1.0,2.0])
print "The theoretical maximum error for the Chebyshev interpolating polynomial of degree 10 is", theoretical_error_chebyshev(10,1.04791*(10**10))

print "The actual maximum error for the cubic splines with n=6 is ", max_error(natural_cubic_spline(f,[-1.0,2.0],6),f,[-1.0,2.0])

print "The actual maximum error for the cubic splines with n=10 is ", max_error(natural_cubic_spline(f,[-1.0,2.0],10),f,[-1.0,2.0])

#Written by Camille Taltas for MT3802

from numpy import *
#Initialize A and B for Exercise 1 
A_1 = matrix([[9,1,2,1,0,4],[0,8,1,2,1,0],[1,0,9,1,2,1],[0,1,0,7,1,2],[0,0,1,0,9,1],[4,0,0,1,0,6]])
B_1 = matrix([[27],[22],[17],[23],[13],[24]])

#Initialize A and B for Exercise 2
A_2 = matrix([[7,2,0,0,0,0],[5,4,3,0,0,0],[0,3,9,1,0,0],[0,0,3,7,1,0],[0,0,0,3,9,1],[0,0,0,0,5,6]])
B_2 = matrix([[20],[34],[47],[27],[17],[17]])

#Initialize funtion that computes omega, the spectral radii and the number of iterations for any 6x6 matrix
def f(A,B, consistently_ordered):
    #Initialize D
    D = zeros((6,6))
    for i in xrange (0,6):
        D[i,i]=A[i,i]

    #Initalize U
    U = zeros((6,6))
    for i in xrange (0,6):
        for j in xrange (0,6):
            if j>i:
                U[i,j]=-A[i,j]

    #Initialize L
    L = zeros((6,6))
    for i in xrange (0,6):
        for j in xrange(0,6):
            if j<i:
                L[i,j] = -A[i,j]

    #Find omega
    if consistently_ordered == True:
        #Finding w, assuming that B_SOR is consistently ordered, 
        B_J = matmul(linalg.inv(D),L+U)
        spectral_rad_J = max(abs(linalg.eigvals(B_J)))
        w = 2.0/(1.0+sqrt(1.0-(spectral_rad_J**2)))
        print "The value of omega is", w
    else :
        #Finding w, with B_SOR not consistently ordered
        spectral_radii = []
        w_vect = []
        w = 1.0 
        while w<2.0:
            B_SOR = matmul(linalg.inv(D-w*L),((1.0-w)*D+w*U))
            spectral_rad = max(abs(linalg.eigvals(B_SOR)))
            spectral_radii.append(spectral_rad)
            w_vect.append(w)
            w = w + 0.01
        w = w_vect[spectral_radii.index(min(spectral_radii))]
        print "The value of omega is", w

    #Spectral Radius of Jacobi Method 
    B_J = matmul(linalg.inv(D),L+U)
    spectral_rad_J = max(abs(linalg.eigvals(B_J)))
    print "The spectral radius of the Jacobi method is", spectral_rad_J

    #Spectral Radius Gauss-Seidel Method
    B_GS = matmul(linalg.inv(D-L),U)
    spec_rad_GS = max(abs(linalg.eigvals(B_GS)))
    print "The spectral radius of the Gauss_Sidel method is", spec_rad_GS

    #Spectral Radius of Successive Relaxaion Method
    B_SOR = matmul(linalg.inv(D-w*L),((1.0-w)*D+w*U))
    spec_rad_SOR = max(abs(linalg.eigvals(B_SOR)))
    print "The spectral radius of the Successive Relaxation method is", spec_rad_SOR

    #Calculating the number of iterations for the Jacobi Method
    x_k = zeros((6,1))
    x_k1 = matmul(B_J,x_k) + matmul(linalg.inv(D),B)
    count = 0 
    while linalg.norm(x_k1-x_k,inf)>=10**(-8):
        x_k = x_k1
        x_k1 = matmul(B_J,x_k) + matmul(linalg.inv(D),B)
        count += 1 
    print "The number of iterations for the Jacobi Method is", count

    #Calculating the number of iterations for the Gauss-Seidel Method
    x_k = zeros((6,1))
    x_k1 = matmul(B_GS,x_k) + matmul(linalg.inv(D-L),B)
    count = 0 
    while linalg.norm(x_k1-x_k,inf)>=10**(-8):
        x_k = x_k1
        x_k1 = matmul(B_GS,x_k) + matmul(linalg.inv(D-L),B)
        count += 1 
    print "The number of iterations for the Gauss_Sidel method is", count

    #Calculating the number of iterations for the Successive Method
    x_k = zeros((6,1))
    x_k1 = matmul(B_SOR,x_k) + matmul(linalg.inv(D-w*L),w*B)
    count = 0 
    while linalg.norm(x_k1-x_k,1)>=10**(-8):
        x_k = x_k1
        x_k1 = matmul(B_SOR,x_k) + matmul(linalg.inv(D-w*L),w*B)
        count += 1 
    print "The number of iterations for the Successive Relaxation method is", count

#Running the function for Exercise 1 
f(A_1,B_1, False)
#Runnnig the fuction for Exercise 2 
f(A_2,B_2,True)

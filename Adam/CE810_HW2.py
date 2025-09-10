import numpy as np
import matplotlib.pyplot as plt

# Define Constants 
EA = 5e7 #N
z = 25 #mm
l = 2500 #mm
Ks = 1.35 #N/mm
del_W = -7 #N
TOL = 1e-4 #N 

# First perform incremental solution 
W = 0
w = 0
N = 0
x_points_1 = [0] #array to store w
y_points_1 = [0] #array to store W 

for i in range(13):
    K_x = EA/l * ((z+w)/l)**2 + Ks + N / l #Update Kg
    w = w + del_W / K_x #update w
    N = EA*((z/l)*(w/l) + 1/2 * (w/l)**2) #update N
    W += del_W #increment W
    x_points_1.append(-w)
    y_points_1.append(-W)

# Restart and preform incremental / iterative solution 
W = 0
w = 0
N = 0
x_points = [0]
y_points = [0]

for i in range(13):
    K_x = EA/l * ((z+w)/l)**2 + Ks + N / l #update Kg
    w = w + del_W / K_x #update w
    N = EA*((z/l)*(w/l) + 1/2 * (w/l)**2) #update N
    W += del_W #increment W 

    # check equilibrium 
    F_unbalanced = N * (z+w)/l + Ks * w - W 
    #print(abs(F_unbalanced))

    #update w until the unbalanced force is less than TOL 
    while abs(F_unbalanced) > TOL:
        K_x = EA/l * ((z+w)/l)**2 + Ks + N / l #update Kg
        w = w - F_unbalanced / K_x #update w 
        N = EA*((z/l)*(w/l) + 1/2 * (w/l)**2) #update N
        F_unbalanced = N * (z+w)/l + Ks * w - W #check unbalanced force 
    
    #print(F_unbalanced)
    #print(W)
    #print(N*(z+w)/l+Ks*w)
        
    x_points.append(-w)
    y_points.append(-W)

#initialize the range of w
x = np.linspace(-60, 60, 5001)
#calculate theoretical solution
y_nospring = -x - 3/2 * x**2 - 1/2 * x**3
y_spring = -3/2 * x - 3/2 * x**2 - 1/2 * x**3

y_2 = (EA/(l**3))*(z**2*x + 1.5*z*x**2 + 0.5*x**3) + Ks*x

#plot the solitions 
plt.figure()
plt.plot(-x, -y_2, label = "Exact Solution" , color = "black")
plt.plot(x_points_1, y_points_1, "--",marker = "x", label = "Incremental Solution", color = "blue")
plt.plot(x_points, y_points, "--",marker = "x", label = "Incremental + Iterative Solution", color = "red")
plt.xlim(0, 60)   
plt.ylim(0, 95) 
plt.ylabel("-W (N)")
plt.xlabel("-w (mm)")
plt.legend() 
plt.grid(True)
plt.show(block = False)

plt.figure() 
plt.plot(-x, y_nospring, label = "Ks = 0", color = "red")
plt.plot(-x, y_spring, label = r"$Ks = \frac{EAz^2}{2l^3}$", color = "blue" )
plt.xlim(0, 2.5)   
plt.ylim(-0.3, 1.6) 
plt.ylabel(r"$\frac{-W}{EA}(\frac{L}{z})^3$")
plt.xlabel(r"$\frac{-w}{z}$")
plt.legend() 
plt.grid(True)
plt.show()

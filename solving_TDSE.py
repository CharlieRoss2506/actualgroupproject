#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: charlotteharrison, charlieross, luzsanchezreal, lucacharalambides
"""

# if using conda in terminal: conda install -c conda-forge findiff
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from findiff import FinDiff
from scipy.sparse.linalg import inv
from scipy.sparse import eye, diags
import matplotlib.animation as animation


# Inputting parameters, defining space & time variables
Nx = 250
xmin = -5
xmax = 5
Nt = 100
tmin = 0
tmax = 60

'''
READ BEFORE INPUTTING POTENTIALS
When inputting own potentials use the variable x_array
Numpy packages can be used to form a string as the potential
'''

# Calculating grid, potential, and initial wave function
x_array = np.linspace(xmin, xmax, Nx)
t_array = np.linspace(tmin, tmax, Nt)


'''
Defining a function to be returned as the default potential 
'''
def f(x_array):
    return x_array**2


def V():
    String = str(input("Please input potential\n"))
    command = """def f(x):
          return """ + String

    exec(command, globals())

    try: 
        return f(x_array)
    except: # add exceptions for errors from program execution.
        print("\nV(x) input does not compute.")
        print("Please input a valid mathematical equation in Numpy format terms of x_array.")
        print("\nSome examples:")
        print("... x_array**2")
        print("... np.sin(x_array)")
        print("\nPlease see the documentation for more details")
        print("\nClosing programme...\n")
        exit()

VV = V()


'''
Converting V into a Diagonal matrix and calculating small psi
'''
Vmatrix = diags(VV)
psi = np.exp(-(x_array+2)**2)

'''
Calculating finite difference elements
'''
dt = t_array[1] - t_array[0]
dx = x_array[1] - x_array[0]

'''
Calculating the Hamiltonian matrix
'''
H = -0.5 * FinDiff(0, dx, 2).matrix(x_array.shape) + Vmatrix 
'''
Note that FinDiff above is a way of representing our partial derivatives. 
FinDiff objects behave like operators; they have a tuple of the form (axis, spacing, degree) in the argument list for each partial derivative.
The last argument stands for the degree of the derivative (in our case 2).
'''

'''
Applying boundary conditions to the Hamiltonian
'''
H[0, :] = H[-1, :] = 0
H[0, 0] = H[-1, -1] = 1

'''
Calculating matrix U, the discretized version of the time propagation operator in the Crank-Nicholson scheme
'''
I_plus = eye(Nx) + 1j * dt / 2. * H # eye returns all elements to 0 except diagonal to 1.
I_minus = eye(Nx) - 1j * dt / 2. * H
U = inv(I_minus).dot(I_plus) 

'''
Iterating over each time, appending each calculation of psi to make a list of all the calculations
'''
psi_list = []
for t in t_array:
    psi = U.dot(psi)
    psi[0] = psi[-1] = 0
    psi_list.append(np.abs(psi))


    
    while True:
        data_input = input("\nWould you like to download a csv file of the data? (y/n)")
        if data_input == 'y':

            df = pd.DataFrame({'x': x, 'psi': psi_list})
            df.to_csv('psi_data.csv')
            print("\nData has been downloaded.")
            break
        elif data_input ==  'n':
            print("\nNo download chosen. ")
            break   
        else:
            print("\nInvalid response.")
            print("Please enter... y (yes for download) or n (no download).")

download_data(Nx, xmin, xmax, tmin, tmax, V()==VV)   
    
    
"""
We now have a numerical solution of simulated data (i.e. `psi_list`) saved to our disk.
In order to visualize this, we plot the data in 2D
"""

fig, ax = plt.subplots()

ax.set_xlabel("x [arb units]")
ax.set_ylabel("$|\Psi(x, t)|$", color="C0")

ax_twin = ax.twinx()
ax_twin.plot(x_array, VV, color="C1") # where time increases vertically on the y-axis. Note the probability density moves back and forth
ax_twin.set_ylabel("V(x) [arb units]", color="C1")

line, = ax.plot([], [], color="C0", lw=2)
ax.grid()
xdata, ydata = [], []

def run(psi):
    line.set_data(x_array, np.abs(psi)**2)
    return line,

ax.set_xlim(x_array[0], x_array[-1])
ax.set_ylim(0, 1)

'''
A better, more cohesive way of viewing the data is with an animation
'''
ani = animation.FuncAnimation(fig, run, psi_list, interval=10)
ani.save("particle_in_a_well.mp4", fps=120, dpi=300) # The code will be saved as an mp4 file

plt.show()

'''
Hence our animation shows the evolution of the probability density function with time, inputting simulated data previously created.
'''

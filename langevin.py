import numpy as np
import matplotlib.pyplot as plt
import math

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from matplotlib.figure import Figure
import matplotlib.animation as animation

import tkinter as Tk



# Initial Conditions
m  = 1		  # Mass (kg)
t0 = 0         # Initial Time (s)
dt = 0.1		 # Time Step (s)
tf = 100		 # Final Time (s)

# Temperature
T = 1
kB = 1

t  = np.linspace(0, tf, int(tf/dt) + 1)	# Time Array

u = 10 			 # Initilal Velocity (m/s)
v = 0			 # Velocity at time t
vAnalytical = [] # Velocity Array (Analytical)
vEuler = [u] 	 # Velocity Array (Euler)

x0 = 0			 # Initial Position (m)
x  = 0			 # Position at time t (m)
xAnalytical = [] # Position Array (Analytical)
xEuler = [x0] 	 # Position Array (Euler)

k = 5			 # V = k*x*x
V = k * x * x	 # Potential Funciton (Nm)
gamma = 1       # Friction Coefficient
sigma = math.sqrt(2 * kB * T * gamma / m)

# Randomly generated numbers
np.random.seed()
N = (tf - t0)/dt
x1 = np.random.uniform(0, 1.0, int(N))
x2 = np.random.uniform(0, 1.0, int(N))
x3 = np.random.uniform(0, 1.0, int(N))
x4 = np.random.uniform(0, 1.0, int(N))

# Calculate acceleration a(x)
def a(x):
	return (-2 * k * x)/m

# Calculate C()
def C():
	c = pow(dt,2) * 0.5 * (a(xEuler[-1]) - gamma * vEuler[-1]) + sigma * math.sqrt(pow(dt,3)) * (0.5 * eta + 0.5 * theta / math.sqrt(3))
	return c


## Analytical Solution
#for i in t:
#	x = x0 + u * i + 0.5 * a * i * i
#	v = u + a * i
#	xAnalytical.append(x)
#	vAnalytical.append(v)

# Euler's Method
for i in t:   
    eta = math.sqrt(-2.0 * math.log(x1[int(i/dt) - 1])) * math.cos(2.0 * math.pi * x2[int(i/dt) - 1]);
    theta = math.sqrt(-2.0 * math.log(x3[int(i/dt) -1])) * math.sin(2.0 * math.pi * x4[int(i/dt) - 1]); 
    if i != 0:
        c = C()
        x = xEuler[-1] + vEuler[-1] * dt + c
        xEuler.append(x)
        v = vEuler[-1] + 0.5 * dt * (a(xEuler[-1]) + a(xEuler[-2])) - dt * gamma * vEuler[-1] + sigma * math.sqrt(dt) * eta - gamma * c
        vEuler.append(v)

# print(len(xEuler))
# print(len(xAnalytical))
# print(len(vEuler))
# print(len(vAnalytical))
# print(len(t))

xEulerArray = np.array(xEuler)
#xAnalyticalArray = np.array(xAnalytical)

vEulerArray = np.array(vEuler)
#vAnalyticalArray = np.array(vAnalytical)


plt.figure(1)
plt.subplot(211)
plt.plot(t, xEulerArray, 'bo') #, xAnalyticalArray, t, 'k')
plt.xlabel('Time -->')
plt.ylabel('Displacement -->')

plt.subplot(212)
plt.plot(t, vEulerArray, 'ro') #, vAnalyticalArray, t, 'r')
plt.xlabel('Time -->')
plt.ylabel('Velocitys -->')
plt.suptitle('Plots for Displacement and Velocity v/s Time')
plt.show()


def animate(i):
    if i < len(xEulerArray):
        a.clear()
        a.plot(xEulerArray, xEulerArray * 0)
        a.plot(xEulerArray[i],0,'ro')


root = Tk.Tk()
root.wm_title("Particle Motion in a Straight Line")


fig = Figure(figsize=(5, 4), dpi=170)
a = fig.add_subplot(111)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)


ani = animation.FuncAnimation(fig, animate, interval = 1)
    
Tk.mainloop()

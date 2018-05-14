import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from matplotlib.figure import Figure
import matplotlib.animation as animation

import tkinter as Tk

# Initial Conditions
m  = 1			 # Mass (kg)
t0 = 0           # Initial Time (s)
dt = 0.1		 # Time Step (s)
tf = 100		 # Final Time (s)

t  = np.linspace(0, tf, int(tf/dt) + 1)	# Time Array


u = 190 			 # Initilal Velocity (m/s)
v = 0			 # Velocity at time t
vAnalytical = [] # Velocity Array (Analytical)
vEuler = [u] 	 # Velocity Array (Euler)

x0 = 0			 # Initial Position (m)
x  = 0			 # Position at time t (m)
xAnalytical = [] # Position Array (Analytical)
xEuler = [x0] 	 # Position Array (Euler)

k = 5			 # V = kx
V = k * x	 	      # Potential Funciton (Nm)
F = - k 		      # Force = -dV/dx (N)
a = F / m 		 # Acceleration (m/s^2)



# Analytical Solution
for i in t:
	x = x0 + u * i + 0.5 * a * i * i
	v = u + a * i
	xAnalytical.append(x)
	vAnalytical.append(v)
    

# Euler's Method
for i in t:
	if i != 0:
		x = xEuler[-1] + vEuler[-1] * dt + 0.5 * a * dt * dt
		v = vEuler[-1] + a * dt
		xEuler.append(x)
		vEuler.append(v)

# print(len(xEuler))
# print(len(xAnalytical))
# print(len(vEuler))
# print(vAnalytical)
# print(len(t))

xEulerArray = np.array(xEuler)
xAnalyticalArray = np.array(xAnalytical)

vEulerArray = np.array(vEuler)
vAnalyticalArray = np.array(vAnalytical)


plt.figure(1)
plt.subplot(211)
plt.plot(t, xEulerArray, 'bo', t, xAnalyticalArray, 'k')
plt.xlabel('Time -->')
plt.ylabel('Displacement -->')

plt.subplot(212)
plt.plot(t, vEulerArray, 'ro', t, vAnalyticalArray, 'k')
plt.xlabel('Time -->')
plt.ylabel('Velocitys -->')
plt.suptitle('Plots for Displacement and Velocity v/s Time')
plt.show()


def animate(i):
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


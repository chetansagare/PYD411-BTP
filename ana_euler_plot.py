import numpy as np
import matplotlib.pyplot as plt

# Initial Conditions
m  = 1		 # Mass (kg)
t0 = 0           # Initial Time (s)
dt = 0.1	 # Time Step (s)
tf = 100	 # Final Time (s)

t  = np.linspace(0, tf, int(tf/dt) + 1)	# Time Array


u = 10 		 # Initilal Velocity (m/s)
v = 0		 # Velocity at time t
vAnalytical = [] # Velocity Array (Analytical)
vEuler = [u] 	 # Velocity Array (Euler)

x0 = 0		 # Initial Position (m)
x  = 0		 # Position at time t (m)
xAnalytical = [] # Position Array (Analytical)
xEuler = [x0] 	 # Position Array (Euler)

k = 5		 # V = kx
V = - k * x	 # Potential Funciton (Nm)
F = k 		 # Force = -dV/dx (N)
a = F / m 	 # Acceleration (m/s^2)



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
plt.ylabel('Distance -->')

plt.subplot(212)
plt.plot(t, vEulerArray, 'ro', t, vAnalyticalArray, 'r')
plt.xlabel('Time -->')
plt.ylabel('Velocity -->')
plt.suptitle('Plots for Distance and Velocity v/s Time')
plt.show()

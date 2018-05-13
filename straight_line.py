import numpy as np
import matplotlib.pyplot as plt

# Initial Conditions
m  = 1			 # Mass (kg)
t0 = 0           # Initial Time (s)
dt = 0.1		 # Time Step (s)
tf = 100		 # Final Time (s)

t  = np.linspace(0, t0, int(tf/dt) + 1)	# Time Array


u = 10 			 # Initilal Velocity (m/s)
v = 0			 # Velocity at time t
vAnalytical = [] # Velocity Array (Analytical)
vEuler = [u] 	 # Velocity Array (Euler)

x0 = 0			 # Initial Position (m)
x  = 0			 # Position at time t (m)
xAnalytical = [] # Position Array (Analytical)
xEuler = [x0] 	 # Position Array (Euler)

k = 5			 # V = kx
V = k * x	 	 # Potential Funciton (Nm)
F = - k * x		 # Force = -dV/dx (N)
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

plt.figure(1)
plt.subplot(211)
plt.plot(xEuler, t, 'bo', xAnalytical, t, 'k')

plt.subplot(212)
plt.plot(vEuler, t, 'ro', vAnalytical, t, 'r')
plt.show()




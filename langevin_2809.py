import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import math

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from matplotlib.figure import Figure
import matplotlib.animation as animation
from operator import add

import tkinter as Tk
import sys
import random
from pylab import show,hist,subplot,figure
import datetime
import sys
import os
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties

def plot_graphs(x_list, y_list, title, xlabel, ylabel):
    directory = 'PLOTS/'
    save_file_name = str(title)
    if not os.path.exists(directory):
        os.makedirs(directory)

    fig = plt.figure(figsize=(8.0, 5.0))
    plt.plot(x_list, y_list, label='title', color='r', figure=fig)
    fontP = FontProperties()
    fontP.set_size('xx-small')
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(str(title))
    plt.legend(prop=fontP)
    file_name = str(directory) + str(save_file_name) + '.pdf'
    plt.savefig(str(file_name), dpi=600)
    plt.clf()
    plt.close()

start_time = datetime.datetime.now()

# Globals 

tau_list = [1,5,10,50,100,500,1000];
steps_in_cycle = 1000
mass = 1
realisations = 1100
transients = 1000
gamma = 1

x0 = 0
v0 = 0

# Spring Constant k(t), V(t) = k(t)*x*x
kmax = 5                        # Max spring constant
kmin =  kmax/2                  # Min spring constant

# Temperature Dependent Contants
T_hot = 4
T_cold = 2
kB = 1
act_hot = 2
act_cold = 1

# Exponential distributed random variable
lambda_hot = 4
lambda_cold = 2
kickV_hot = 1
kickV_cold = 0.5


print("All constants compiled.")

# Sigma Function
def sigma(i,tau):
    return math.sqrt(2 * kB * T(i,tau) * gamma / mass)

# Temperature function 
def T(i,tau):
    i = i%tau
    if i >= tau/2:
        return T_hot
    else:
        return T_cold

# Spring constant as function of time 
def k(i,tau):
    i = i%tau
    kslope = (kmax - kmin)/(tau/2)
    if i<tau/2:
        return (kmin - kslope * (i - tau/2))
    else:
        return (kmax + kslope * (i - tau))

# dk/dt as a function of time 
def kdot(i,tau):
    i = i%tau
    kslope = (kmax - kmin)/(tau/2)
    if i >= tau/2:
        return kslope
    else:
        return (-kslope)

# Calculate acceleration a(x, t)
def a(x, i,tau):
    return ( -2 * k(i,tau) * x)/mass

# Calculate C()
def C(i, xlast, vlast, eta, theta,tau,dt):
    return pow(dt,2) * 0.5 * (a(xlast, i,tau) - gamma * vlast) + sigma(i,tau) * math.sqrt(pow(dt,3)) * (0.5 * eta + 0.5 * theta / math.sqrt(3))
    

# Position Function
def X(xlast,vlast,dt,c):
    return xlast + vlast * dt + c

# Velocity Function 
def V(vlast, xlast, dt, x_sec_last, c, eta, tau, t):
    return vlast + 0.5 * dt * ( a(xlast, t+1, tau) + a(x_sec_last, t, tau)) - dt * gamma * vlast + sigma(t, tau) * math.sqrt(dt) * eta - gamma * c

def W(wlast, t,xlast, dt, tau):
    return wlast + (0.5 * kdot(t,tau) * xlast * xlast * dt)

def Q(qlast, vlast, eta, dt):
    return qlast + ( -gamma * vlast * vlast + eta * vlast) * dt

# Lambda function
def lamda(i,tau):
    i = i%tau
    if i >= tau/2:
        return lambda_hot
    else:
        return lambda_cold

# Kick Velocity function
def kickVel(i,tau):
    i = i%tau
    direction = random.uniform(0,1)
    if direction < 0.5:
        direction = -1
    else:
        direction = +1
    if i >= tau/2:
        return (kickV_hot * direction)
    else:
        return (kickV_cold * direction)

w_tau_list = []
q_tau_list = []

for tau in tau_list[5:6]:
    print("Tau : "+str(tau))
    dt = tau / steps_in_cycle
    time_list = np.linspace(0,tau * realisations, steps_in_cycle * realisations +1)
    eta_list = np.random.normal(0, 1.0, len(time_list)+1)
    theta_list = np.random.normal(0, 1.0, len(time_list)+1)

    x_list = []
    v_list = []
    x_sqr_list = []
    v_sqr_list = []
    w_list = []
    q_list = []
    k_list = []
    
    x = x0
    v = v0
    w = 0
    q = 0
    count = 0
    
    for t in time_list:
        if t == 0:
            timeKick_uni = random.uniform(0, 1)
            timeKick = (-1 / lamda(t, tau)) * math.log(abs(1 - timeKick_uni)) + t
            continue
        eta = eta_list[int(t/dt)]
        theta = theta_list[int(t/dt)]
        c = C(t,x,v,eta,theta,tau,dt)
        xlast = x
        vlast = v
        x = X(x, v, dt, c)
        v = V(v, x, dt, xlast, c, eta, tau, t)
        wlast = w
        qlast = q
        w = W(w, t, x,dt,tau)
        q = Q(q, v, eta, dt)
        
        # kicks
        if t >= timeKick:
            v = v + kickVel(t, tau)
            
            timeKick_uni = random.uniform(0, 1)
            timeKick = (-1 / lamda(t, tau)) * math.log(abs(1 - timeKick_uni)) + t
            
            

        if t / tau > transients:
            
            if len(x_list)==0:
                x_list.append(x)
                x_sqr_list.append(x*x)
                v_list.append(v)
                v_sqr_list.append(v*v)
                count += 1
            else:
                x_list.append((x + x_list[-1] * count) / (count + 1))
                x_sqr_list.append((x*x + x_sqr_list[-1] * count)/ (count + 1))
                v_list.append((v + v_list[-1] * count)/ (count + 1))
                v_sqr_list.append((v*v + v_sqr_list[-1] * count)/ (count + 1))
                count += 1

            k_list.append(k(t,tau))
            w_list.append(w)
            q_list.append(q)

    slice_list = np.linspace(0,steps_in_cycle,len(x_list[-steps_in_cycle:]))
    
    plot_graphs(slice_list, x_list[-steps_in_cycle:],'x_vs_t','t','x')
    plot_graphs(slice_list, x_sqr_list[-steps_in_cycle:],'x_sqr_vs_t','t','x2')
    plot_graphs(slice_list, v_list[-steps_in_cycle:],'v_vs_t','t','v')
    plot_graphs(slice_list, v_sqr_list[-steps_in_cycle:],'v_sqr_vs_t','t','v2')
    plot_graphs(slice_list, w_list[-steps_in_cycle:],'w_vs_t','t','w')
    plot_graphs(slice_list, q_list[-steps_in_cycle:],'q_vs_t','t','q')
    plot_graphs(slice_list, k_list[-steps_in_cycle:],'k_vs_t','t','k')

    w_tau_list.append(np.mean(w_list))
    q_tau_list.append(np.mean(q_list))
    
print("Loops completed started mean")


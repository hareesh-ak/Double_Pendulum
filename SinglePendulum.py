# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 13:01:16 2021

@author: Hareesh
"""


import numpy as np
import matplotlib.pyplot as plt
import time 
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

g=10

class Pendulum(object):
    mass= 1
    def __init__(self,length,theta):
        #self.mass = mass
        #self.length = length
        self.theta=theta*np.pi/180
        self.x_init = length*np.sin(self.theta)
        self.y_init = -length*np.cos(self.theta)
        self.length= length
    thetadot=0
    fig=plt.figure(1,figsize=(11,10))
    axis=plt.axes()
    def draw_fig(self):
        Axis=self.axis
        w=10;h=10
        Axis.set_ylim(-h-self.length,h+self.length)
        Axis.set_xlim(-w-self.length,w+self.length)
        Axis.scatter(self.x_init,self.y_init)
        #m=self.initial_pos[1]/self.initial_pos[0]
        Axis.plot([0,self.x_init],[0,self.y_init])
        time.sleep(0.001)
        
        
    def update_pos(self,dt):
        z=update_step([self.theta,self.thetadot],dt)
        self.theta,self.thetadot = z[0],z[1]
        self.x_init = self.length*np.sin(self.theta)
        self.y_init = -self.length*np.cos(self.theta)
                
        
###
### Input length and inital angle of pendulum below:
### p is an instance of the Pendulum class
###
p=Pendulum(length=10,theta=-90) #theta in degress from the vertical axis
p.draw_fig()

##
##  update_step used to update theta and theta_dot of pendulum; 
##  fun is used in the ODE solver solve_ivp from scipy.integrate
##


def fun(t,X):   
    x,y = X
    xdot = y
    ydot = -(g/p.length)*np.sin(x)
    return xdot, ydot

def update_step(init_pos,dt):
    X0=np.array([init_pos[0],init_pos[1]])
    z=solve_ivp(fun,[0,dt],X0,method='RK45',first_step=dt)
    return z.y[0][-1],z.y[1][-1]


#p.update_pos(.1)

# h s the time step 
# animate is defined for FuncANimation from matplotlib.animation    

h=.1
def animate(i):
   p.axis.cla()
   p.update_pos(h)
   p.draw_fig()
   
   

fig=plt.figure(1)
anime = FuncAnimation(fig,animate,interval=100)
plt.show()

    
    
   
    



    
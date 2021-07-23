# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 20:12:55 2021

@author: Hareesh
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

plt.style.use('dark_background')
pi=np.pi
s=time.time()
g=9.81
"""
The functions defined below are called while evolving the system of differential equations
used to decribe the dynamics of double pendulum. Check report for details about these functions

###INSERT REPORT LINK HERE!!
"""



    

class Double_Pendulum(object):
    m1=1;m2=1
    w1=0;w2=0
    t=0
    def __init__(self,l1=1,l2=1,th1=180,th2=180.1,color='b'):
        """
        

        Parameters
        ----------
        l1 : float
            length of string attached to mass m1.
        l2 : float
            length of string attached to mass m2.
        th1 : float
            angle(in degrees) made by string l1 with vertical.
        th2 : float
            angle(in degrees) made by string l2 with vertical.
        color : string
            The supported color abbreviations are the single letter codes;
            passed as arg for matplotlib.pyplot.plot( .scatter)
    
            =============    ===============================
            character        color
            =============    ===============================
            ``'b'``          blue
            ``'g'``          green
            ``'r'``          red
            ``'c'``          cyan
            ``'m'``          magenta
            ``'y'``          yellow
            ``'k'``          black
            ``'w'``          white
            =============    ===============================

        """
        #self.mass = mass
        #self.length = length
        self.th1=th1*np.pi/180
        self.th2=th2*np.pi/180
        self.x1_init = l1*np.sin(self.th1)
        self.y1_init = -l1*np.cos(self.th1)
        self.x2_init = l1*np.sin(self.th1) + l2*np.sin(self.th2)
        self.y2_init = -(l1*np.cos(self.th1) + l2*np.cos(self.th2))
        self.l1= l1;self.l2= l2
        self.color=color
        self.init_pos=np.array([self.th1,self.th2,self.w1,self.w2])
    
    fig,axis=plt.subplots()#,figsize=(11,10))
    axis=plt.axes()
    
    
    def draw_fig(self):
        """
        

        Returns
        -------
        return Nonetype.
        Plots current configuration of system in matplotlib figure.

        """
        Axis=self.axis
        w=1;h=1
        Axis.set_ylim(-h-self.l1-self.l2,h+self.l1+self.l2)
        Axis.set_xlim(-w-self.l1-self.l2,w+self.l1+self.l2)
        Axis.scatter(self.x1_init,self.y1_init,c=self.color)
        Axis.scatter(self.x2_init,self.y2_init,c=self.color)
        #m=self.initial_pos[1]/self.initial_pos[0]
        Axis.plot([0,self.x1_init],[0,self.y1_init],c=self.color)
        Axis.plot([self.x1_init,self.x2_init],[self.y1_init,self.y2_init],self.color)
        Axis.set_xlabel(f'Time = {round(self.t,2)} sec',fontsize=15)
        #time.sleep(0.001)
        
        
    def A1(self,th1,th2,w1,w2):
        return (self.m2*self.l2)/((self.m1+self.m2)*self.l1)*np.cos(th1-th2)

    def A2(self,th1,th2,w1,w2):
        return self.l1/self.l2*np.cos(th1-th2)
    
    def B1(self,th1,th2,w1,w2):
        return (self.m2*self.l2)/self.l1/(self.m1+self.m2)*w2**2*np.sin(th1-th2) + \
            g/self.l1*np.sin(th1)
    
    def B2(self,th1,th2,w1,w2):
        return -self.l1/self.l2*w1**2*np.sin(th1-th2) + g/self.l2*np.sin(th2)
        
    
    def F1(self,th1,th2,w1,w2):
        return (-self.B1(th1,th2,w1,w2) +\
                self.A1(th1,th2,w1,w2)*self.B2(th1,th2,w1,w2))/ \
            (1-self.A1(th1,th2,w1,w2)*self.A2(th1,th2,w1,w2))
    
    def F2(self,th1,th2,w1,w2):
        return (-self.B2(th1,th2,w1,w2) +self.A2(th1,th2,w1,w2)* \
                self.B1(th1,th2,w1,w2))/(1-self.A1(th1,th2,w1,w2)* \
                                                        self.A2(th1,th2,w1,w2))
    
    
    def f(self,t,X):
        th1,th2,w1,w2=X
        th1dot = w1;th2dot=w2
        w1dot = self.F1(th1,th2,w1,w2)
        w2dot = self.F2(th1,th2,w1,w2)
        return th1dot,th2dot,w1dot,w2dot
    
    locations=np.array([[],[],[],[]])
    """
    The function used to evolve system of equations for a single time step.
    """
    def update_step(self,dt):
        X0=self.init_pos
        z=solve_ivp(self.f,[0,dt],X0,method='RK45',first_step=dt)
        self.init_pos=z.y[0][-1],z.y[1][-1],z.y[2][-1],z.y[3][-1]
        self.t+=dt
        #self.locations=np.append(self.locations,[[self.init_pos[0]],[self.init_pos[1]],[self.init_pos[2]],[self.init_pos[3]]])
        return self.init_pos
    
    def get_state(self):
        print(self.th1,self.th2,self.w1,self.w2)
    
    def update_pos(self,dt):
        """
        

        Parameters
        ----------
        dt : float
            time step used for evolving system;
            passed as arg for scipy.integrate.solve_ivp

        Returns
        -------
        None.

        """
        z=self.update_step(dt)
        self.th1,self.th2,self.w1,self.w2 = z
        self.x1_init = self.l1*np.sin(self.th1)
        self.y1_init = -self.l1*np.cos(self.th1)
        self.x2_init = self.l1*np.sin(self.th1) + self.l2*np.sin(self.th2)
        self.y2_init = -(self.l1*np.cos(self.th1) + self.l2*np.cos(self.th2))
    def get_trajectory(self,dt,tf):
        X0=self.init_pos
        z=solve_ivp(self.f,t_span=(0,tf),y0=X0,first_step=dt,max_step=dt)#,min_step=dt)
        
        
        return z.t,z.y

p1=Double_Pendulum(l1=1,l2=1,th1=90.,th2=90.2)
p2=Double_Pendulum(l1=1,l2=1,th1=90.,th2=90.3,color='r')

#p1.draw_fig()
#p2.draw_fig()

h=.1
def animate(i):
   p1.axis.cla()
   p2.axis.cla()   
   p1.update_pos(h)
   p2.update_pos(h)
   p1.draw_fig()
   p2.draw_fig()
   
   

fig=plt.figure(1)
plt.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)
anime = FuncAnimation(fig,animate,frames=300,interval=100,cache_frame_data=False)
#plt.close()
anime.save('trial.mp4',writer='ffmpeg',bitrate=400,dpi=300)

"""
plt.show()
plt.style.use('default')

#init1=np.array([0,0,0,0])
#init2=np.array([0,0,0,0])

t1,soln1=p1.get_trajectory(dt=0.01,tf=20*pi) 
#t2,soln2=p2.get_trajectory(dt=0.001,tf=10*pi)
e=time.time()

print((e-s)/60.0)
"""         
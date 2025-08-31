# -*- coding: utf-8 -*-
"""
CE 810, Homework 2
Created on Sat Aug 30 22:56:16 2025

@author: Xinlong Du, University of Kansas
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import numpy as np

# set the font parameters
mpl.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 1
tick_width = 1
tick_mj_sz = 2
tick_mn_sz = 1
plt_line_width = 0.8

small_fig_size = (3.5,3.5)
big_fig_size = (6,5)
fig_font_size = 8

def force_disp_plot(x,y,file_name):

    fig = plt.figure(figsize=small_fig_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_tick_params(which='major', size=tick_mj_sz, width=tick_width, direction='in')
    ax.xaxis.set_tick_params(which='minor', size=tick_mn_sz, width=tick_width, direction='in')
    ax.yaxis.set_tick_params(which='major', size=tick_mj_sz, width=tick_width, direction='in')
    ax.yaxis.set_tick_params(which='minor', size=tick_mn_sz, width=tick_width, direction='in')
    
    # plot
    ax.plot(x,y, linewidth=plt_line_width)
    plt.grid()
    
    # Add the x and y-axis labels
    ax.set_xlabel('-w/z')
    ax.set_ylabel('-Wl^3/EAz^3')
    
    # Save figure
    plt.savefig('./'+file_name+'.svg', transparent=False, bbox_inches='tight')
    
def force_disp_plot2(EA,z,L,Ks,x,y,file_name):

    fig = plt.figure(figsize=small_fig_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_tick_params(which='major', size=tick_mj_sz, width=tick_width, direction='in')
    ax.xaxis.set_tick_params(which='minor', size=tick_mn_sz, width=tick_width, direction='in')
    ax.yaxis.set_tick_params(which='major', size=tick_mj_sz, width=tick_width, direction='in')
    ax.yaxis.set_tick_params(which='minor', size=tick_mn_sz, width=tick_width, direction='in')
    
    # plot
    ax.plot(x,y, linewidth=plt_line_width)
    x2=-np.linspace(0.0,60,300);
    y2=EA/L**3*(z**2*x2+1.5*z*x2**2+0.5*x2**3)+Ks*x2;
    ax.plot(-x2,-y2, linewidth=plt_line_width)
    plt.grid()
    
    # Add the x and y-axis labels
    ax.set_xlabel('-w/z')
    ax.set_ylabel('-Wl^3/EAz^3')
    ax.legend(['Incremental', 'Exact'])
    
    # Save figure
    plt.savefig('./'+file_name+'.svg', transparent=False, bbox_inches='tight')

#%% Problem 1    
x=np.linspace(0.0,2.4,100);
y=x-1.5*x**2+0.5*x**3;
force_disp_plot(x,y,'HW2P1')

#%% Problem 2
x=np.linspace(0.0,2.1,100);
y=1.5*x-1.5*x**2+0.5*x**3;
force_disp_plot(x,y,'HW2P2')

#%% Problem 3
EA=5e7;
z=25;
L=2500;
Ks=1.35;
dW=-7;

#initial condition
w=0;
N=0;
wList=[];
WList=[];
for i in range(0,14):
    wList.append(w);
    WList.append(dW*i);
    Kt=EA/L*(z+w)**2/L**2+N/L+Ks;
    w=w+dW/Kt;
    N=EA*(z/L*w/L+0.5*(w/L)**2);
force_disp_plot2(EA,z,L,Ks,-np.array(wList),-np.array(WList),'HW2P3')

#%% Problem 4
tol=1e-4; #tolerance for Newton-Rapson iteration
#initial condition
w=0;
N=0;
wList=[];
WList=[];
for i in range(0,14):
    wList.append(w);
    WList.append(dW*i);
    Kt=EA/L*(z+w)**2/L**2+N/L+Ks;
    w=w+dW/Kt;
    N=EA*(z/L*w/L+0.5*(w/L)**2);
    g=N*(z+w)/L+Ks*w-dW*(i+1);
    while abs(g)>tol:
        Kt=EA/L*(z+w)**2/L**2+N/L+Ks;
        w=w-g/Kt;
        N=EA*(z/L*w/L+0.5*(w/L)**2);
        g=N*(z+w)/L+Ks*w-dW*(i+1);
force_disp_plot2(EA,z,L,Ks,-np.array(wList),-np.array(WList),'HW2P3')

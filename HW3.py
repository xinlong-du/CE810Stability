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
    x2=-np.linspace(0.0,60,300);
    y2=EA/L**3*(z**2*x2+1.5*z*x2**2+0.5*x2**3)+Ks*x2;
    ax.plot(-x2,-y2, linewidth=plt_line_width)
    ax.plot(x,y, linewidth=plt_line_width)
    plt.grid()
    
    # Add the x and y-axis labels
    ax.set_xlabel('-w')
    ax.set_ylabel('-W')
    ax.legend(['Exact','Numerical'])
    
    # Save figure
    plt.savefig('./'+file_name+'.svg', transparent=False, bbox_inches='tight')

#%% Problem 4
EA=5e7;
z=25;
L=2500;
Ks=1.35;
tol=1e-4; #tolerance for Newton-Rapson iteration

#initial condition
p=np.array([0,0]); #[u,w]
N=0;

#force increment vector dW=-7; dU=70;
dq=np.array([70,-7])

wList=[];
WList=[];
for i in range(0,17):
    wList.append(p[1]);
    WList.append(dq[1]*i);
    
    beta=(z+p[1])/L;
    Kt=EA/L*np.array([[1,-beta],[-beta,beta**2+Ks*L/EA]])+np.array([[0,0],[0,N/L]]);
    p=p+np.matmul(np.linalg.inv(Kt),dq);
    N=EA*(-p[0]/L+z/L*p[1]/L+0.5*(p[1]/L)**2);
    g=N/L*np.array([-L,z+p[1]])+np.array([0,Ks*p[1]])-dq*(i+1);
    while np.linalg.norm(g)>tol:
        Kt=EA/L*np.array([[1,-beta],[-beta,beta**2+Ks*L/EA]])+np.array([[0,0],[0,N/L]]);
        p=p-np.matmul(np.linalg.inv(Kt),g);
        N=EA*(-p[0]/L+z/L*p[1]/L+0.5*(p[1]/L)**2);
        g=N/L*np.array([-L,z+p[1]])+np.array([0,Ks*p[1]])-dq*(i+1);
force_disp_plot2(EA,z,L,Ks,-np.array(wList),-np.array(WList),'HW3P1')

# -*- coding: utf-8 -*-
"""
CE 810, Homework 3
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
    
def force_disp_plot2(x,y,hon,ver,file_name):

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
    ax.set_xlabel(hon)
    ax.set_ylabel(ver)
    
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

uList=[];
UList=[];
wList=[];
WList=[];
for i in range(0,38):
    uList.append(p[0]);
    UList.append(dq[0]*i);
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
        beta=(z+p[1])/L;
        N=EA*(-p[0]/L+z/L*p[1]/L+0.5*(p[1]/L)**2);
        g=N/L*np.array([-L,z+p[1]])+np.array([0,Ks*p[1]])-dq*(i+1);
force_disp_plot2(np.array(uList),np.array(UList),'u','U','HW3P1')
force_disp_plot2(-np.array(wList),-np.array(WList),'-w','-W','HW3P2')

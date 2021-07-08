# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:08:18 2020

@author: Charlotte
"""

from sympy.solvers import solve, solveset
from sympy import Symbol, sin, cos, pi, nsolve
import numpy as np
from matplotlib import pyplot as plt

v1m = Symbol('v1m')
v1s = Symbol('v1s')
theta_m = Symbol('theta_m')
pars = [v1m, v1s, theta_m]
x0 = [1,1,1]

theta_i = pi
mm = 4
ms = 195
v0m = 1
#theta_s = 0.25*pi+pi



def main():   
    solutionslist=[]
    for theta_s in np.arange(pi, 1.250001*pi, 0.01*pi):
        print(theta_s)
        eqs = [sin(theta_m)*mm*v1m+sin(theta_s)*ms*v1s-sin(theta_i)*mm*v0m,
               cos(theta_m)*mm*v1m+cos(theta_s)*ms*v1s-cos(theta_i)*mm*v0m,
               mm*v1m**2+ms*v1s**2-mm*v0m**2]
        
        solutions = solve(eqs, pars)
#        print ("solutions",solutions,'\n')
    
        solutions_filtered=filter_trivial(solutions, theta_s)
#        print("solutions_filtered", solutions_filtered)
    
        sol = find_solution(solutions_filtered, theta_s)
        solutionslist.append((*sol, theta_s))
#    print ('final solution', solutionslist)
        
    np.savetxt('solutionslist.txt', np.array([*solutionslist]),header='v1m/v0m, v1s/v0m, theta_m, theta_s')
    
    plot_collisions(solutionslist, save=True)

    
def filter_trivial(solutions, theta_s):
    """
    Filters out:
    v1m = v0m (either because v1m=v0m or v1m=-v0m and angle is flipped)
    v2m negative
    v1m = 0 and theta_s!=pi (v1m can be 0 but only if theta_s=theta_i)
    """
    solutions_filtered=[]
    for sol in solutions:
        if np.absolute(sol[0]) == v0m:
            continue
        if sol[1] < 0:
            continue
        if (sol[0]==0 and theta_s!=theta_i):
            continue      
        solutions_filtered.append(sol)    
    return solutions_filtered

def find_solution(solutions, theta_s):
    """
    Assumes that there are no trivial solutions in the list
    if v1m is 0, all solutions for theta_m are ok. but this only works if theta_s=theta_i. (but that is checked in filter_trivial)

    sin(theta_m) should have opposite sign of sin(theta_s)  (or both zero) - ONLY FOR THETA_I=PI 
    v1m and v1s should be positive (or zero)
    
    if there is no solution, convert wrong solution into right one by:
    changing sign of v1m
    flipping angle (adding pi)  

    """
    for sol in solutions:
        if sol[0]==0:
#            print(sol)
            sol = (sol[0], sol[1], 0)#set angle to 0 by default
            return sol
            
        if (((sin(sol[2])==sin(theta_s)==0) or (sin(sol[2])!=0 and sin(theta_s)!=0 
              and sin(sol[2])*sin(theta_s)<0)) and sol[0]>0 and sol[1]>=0):
#            print(sol)
            return sol
    return (-sol[0], sol[1], sol[2]-pi)
            

    


def plot_arrows(vm, vs, theta_m, theta_s, fig, ax):
    r=1
    base_x = float(r*sin(theta_s-pi))
    base_y = float(np.absolute(r*cos(theta_s-pi)))
#    print(base_x, base_y)
    
    m_x = float(vm/v0m*sin(theta_m))
    m_y = float(vm/v0m*cos(theta_m))
    
    s_x = float(vs/v0m*sin(theta_s))
    s_y = float(vs/v0m*cos(theta_s))
    

    if vm!=0:
        ax.arrow(base_x, base_y, m_x, m_y, length_includes_head=True, width=0.01)
    else:
        print("vm =0!")
    if vs!=0:
        ax.arrow(base_x, base_y, s_x, s_y, length_includes_head=True, width=0.01)
    else:
        print("vs =0!")
    

def plot_collisions(solutionslist, save=False):
    fig, ax = plt.subplots()
    ax.axis('equal')
    ax.axis([-0.5,2,-0.5,2])
    circ1 = plt.Circle((0,0),1)
    ax.add_artist(circ1)
    for sol in solutionslist:
        plot_arrows(*sol, fig, ax)
    ax.set_aspect('equal', 'box')
    fig.savefig('collisions.png')
    fig.show()
    

main()

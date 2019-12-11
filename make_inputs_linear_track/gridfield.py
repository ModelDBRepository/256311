#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 09:03:49 2018

@author: spiros
"""

def gridfield(theta, lambda_var, xo, yo, x, y):
    import numpy as np
    '''
    Description goes here
    '''
    th1 = np.array([np.cos(theta), np.sin(theta)]).reshape(-1,1)
    th2 = np.array([np.cos(theta + np.pi/3), np.sin(theta + np.pi/3)]).reshape(-1,1)
    th3 = np.array([np.cos(theta + 2*np.pi/3), np.sin(theta + 2*np.pi/3)]).reshape(-1,1)
    
    x -= xo
    y -= yo
    
    y /= float(lambda_var)
    x /=float(lambda_var)
    
    p = np.array([x,y]).reshape(-1,1)
    
    g = (1/4.5) * ( np.cos(np.dot(p.T, th1)) + np.cos(np.dot(p.T, th2)) + np.cos(np.dot(p.T, th3)) + 1.5).item()

    return g
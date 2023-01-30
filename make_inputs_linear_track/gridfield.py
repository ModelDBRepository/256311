#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 09:03:49 2018.

@author: spiros
"""
import numpy as np


def gridfield(theta, lambda_var, xo, yo, x, y):
    """
    Generate grid-like fields in a specific location.

    Parameters
    ----------
    theta : FLOAT
        DESCRIPTION.
    lambda_var : FLOAT
        DESCRIPTION.
    xo : INTEGER
        x-coordinate of the place field.
    yo : INTEGER
        y-coordinate of the place field.
    x : INTEGER
        x point in space.
    y : INTEGER
        y point in space.

    Returns
    -------
    g : FLOAT
        firing rate at point x,y based on place field (xo, yo).

    """
    th1 = np.array([np.cos(theta), np.sin(theta)]).reshape(-1, 1)
    th2 = np.array(
        [np.cos(theta + np.pi/3), np.sin(theta + np.pi/3)]).reshape(-1, 1)
    th3 = np.array([np.cos(theta + 2*np.pi/3),
                    np.sin(theta + 2*np.pi/3)]).reshape(-1, 1)

    x -= xo
    y -= yo

    y /= float(lambda_var)
    x /= float(lambda_var)

    p = np.array([x, y]).reshape(-1, 1)

    g = (1/4.5) * (np.cos(np.dot(p.T, th1)) +
                   np.cos(np.dot(p.T, th2)) + np.cos(np.dot(p.T, th3)) +
                   1.5).item()

    return g

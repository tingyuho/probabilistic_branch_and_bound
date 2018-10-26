# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 22:35:17 2017

@author: TingYu Ho
"""
import math
import numpy as np

def testfunction_shifted_sin (X,r):
    mu, sigma = 0, 0 # mean and standard deviation
    s = np.random.normal(mu, sigma, r)
    x1 = X[0]
    x2 = X[1]
    noise_funtion = -2.5*math.sin(math.pi*float(x1+60)/180)*math.sin(math.pi*float(x2+60)/180)-math.sin(math.pi*float(x1+60)/(36))*math.sin(math.pi*float(x2+60)/(36)) + s
    l_noise_funtion = noise_funtion.tolist()
    return l_noise_funtion

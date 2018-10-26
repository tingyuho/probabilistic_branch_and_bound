# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:47:31 2017

@author: TingYu Ho
r: number of replication
"""
import numpy as np

def testfunction_rosenbrock(X,r):
    mu, sigma = 0, 0 # mean and standard deviation
    s = np.random.normal(mu, sigma, r)
    x = X[0]
    y = X[1]
    a = 1. - x
    b = y - x*x
    noise_funtion=a*a + b*b*100
    l_noise_funtion = noise_funtion.tolist()
    return l_noise_funtion

if __name__ == "__main__":
    print(testfunction_rosenbrock([1,2], 2))
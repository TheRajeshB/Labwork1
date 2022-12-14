# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:54:24 2021

@author: farzh
"""
import numpy as np

#calculating standard deviation in 2-passes

#loading data into array 
data_1 = np.loadtxt('cdata.txt')


sum = 0.0 

#loop to calculate sum of values
for i in data_1:
    sum += i

#calculating mean    
mean = sum/len(data_1)
print('mean to be used in eqn 1 is:',mean)

#calculating standard deviation
sum_diff = 0.0

#loop to calculate (x_i - mean(x))**2 and find the sum
for j in data_1:
    diff = (j - mean)**2 #square of difference between data points and mean of data points
    sum_diff = sum_diff + diff

#standard deviation cal
sigma_2pass = (sum_diff/(len(data_1)-1))**0.5
print('standard deviation using eqn 1 is:',sigma_2pass)

#calculating standard deviation using 1-pass of data 

sum_2 = 0.0 
sum_squared = 0.0
n = len(data_1)

for k in data_1:
    sum_2 += k #sum of all data points
    x_i_squared = k**2 #calulating square of each data point
    sum_squared += x_i_squared #calculating sum of square of all data points
    
    
   



mean_2 = sum_2/len(data_1) #mean to be used in eqn 2
print('mean to be used in eqn2 is:', mean_2)

#calculating standard deviation

#check to see if a negative sign is incurred under the square root for eqn2
if sum_squared - n*(mean_2)**2  < 0:
    print('there is a negative sign under the square root for eqn 2')
else: 
    sigma_1pass = ((sum_squared - n*(mean_2)**2)/(n-1))**0.5
print('standard deviation using eqn 2 is:',sigma_1pass)

print('Standard deviation from numpy method is:',np.std(data_1, ddof=1))
#calculating relative error of values from eqns 1 and 2 with respect to numpy method 

relative_error_eqn1 = (sigma_2pass - np.std(data_1, ddof=1))/(np.std(data_1, ddof=1))
relative_error_eqn2 = (sigma_1pass - np.std(data_1, ddof=1))/(np.std(data_1, ddof=1))
print('Relative error for eqn 1 is:',relative_error_eqn1)
print('Relative error for eqn2 is:', relative_error_eqn2)

if np.abs(relative_error_eqn1) < np.abs(relative_error_eqn2):
    print('eqn 2 has larger relative error by a factor of:',np.abs(relative_error_eqn2/relative_error_eqn1 ) )
else:
    print('eqn 1 has larger relative error by a factor of:',np.abs( relative_error_eqn1/relative_error_eqn2))

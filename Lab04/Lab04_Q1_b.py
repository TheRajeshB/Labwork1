# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 22:59:37 2021

@author: farzh
"""
from numpy.random import rand
from numpy.linalg import solve
from numpy import dot,mean
from Lab04_functions import *
from time import time
import matplotlib.pyplot as plt

#lists to store values
 
err_Gauss_list = [] #list to store error for Gaussian elimination method
err_PP_list = [] #list to store error for Partial Pivoting method
err_LU_list =[]#list to store error for LU method

Gauss_time=[] #list to time take to find solution using Gaussian elimination method
PP_time = [] #list to time take to find solution using Partial Pivot method
LU_time = [] #list to time take to find solution using LU decomposition method
N_values = []

for N in range(5,300): #loop to form arrays and vectors for various N values
    N_values.append(N)
    v = rand(N)
    A = rand(N, N)
    
    start_Gauss = time()
    x_Gauss = GaussElim(A,v) #Gaussian elimination method
    end_Gauss = time()
    
    start_PP = time()
    x_PP = PartialPivot(A,v) #Partial pivoting method
    end_PP = time()
   
    start_LU = time()
    x_LU = solve(A, v) #LU decomposition method 
    end_LU = time() 

    
    v_sol_Gauss = dot(A,x_Gauss) #results of product of A and x, with x using Gaussian method
    v_sol_PP= dot(A,x_PP)        #results of product of A and x, with x using Partial Pivot method 
    v_sol_LU = dot(A, x_LU)      #results of product of A and x, with x using LU decomposition method
    
    
    err_Gauss = mean(abs(v-v_sol_Gauss)) # error using Gaussian elimination method
    err_PP    = mean(abs(v-v_sol_PP))    # error using Partial Pivoting method
    err_LU    = mean(abs(v-v_sol_LU))    # error using LU decomposition method
    
    
    #adding values to corresponding lists
    err_Gauss_list.append(err_Gauss)
    err_PP_list.append(err_PP)
    err_LU_list.append(err_LU)
    
    
    
    Gauss_time.append(end_Gauss- start_Gauss)
    PP_time.append(end_PP - start_PP)
    LU_time.append(end_LU - start_LU)


plt.errorbar(N_values,err_Gauss_list,label='error using Gaussian elimination')
plt.errorbar(N_values,err_PP_list,label='error using Partial Pivoting')
plt.errorbar(N_values,err_LU_list,label='error using LU decomposition')
plt.title('Comparison of errors using Gaussian Elim, Partial Pivoting and LU decomp')
plt.yscale('log')
plt.xlabel('N values')
plt.ylabel('error onn log scale')
plt.legend()
plt.show()

plt.scatter(N_values,Gauss_time,marker='o',label='Time taken by Gaussian elimination method')
plt.scatter(N_values,err_PP_list,marker='o',label='Time taken by Partial Pivoting method')
plt.scatter(N_values,LU_time,marker='o',label='Time taken by LU decomposition method')
plt.title('Comparison of time taken by Gaussian Elim, Partial Pivoting and LU decomp')
plt.yscale('log')
plt.xlabel('N values')
plt.ylabel('Time taken on log scale')
plt.legend()
plt.show()

plt.errorbar(N_values,Gauss_time,label='Time taken by Gaussian elimination method')
plt.errorbar(N_values,err_PP_list,label='Time taken by Partial Pivoting method')
plt.errorbar(N_values,LU_time,label='Time taken by LU decomposition method')
plt.title('Comparison of time taken by Gaussian Elim, Partial Pivoting and LU decomp')
plt.yscale('log')
plt.xlabel('N values')
plt.ylabel('Time taken on log scale')
plt.legend()
plt.show()
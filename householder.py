import sys
import os
import numpy as np

from reflection import Reflection

class Householder:
    """Performs QR decomposition of a matrix A to QR using the Housholder transformation.
    
    ...
    
    
    Attributes
    ----------
    A: numpy matrix
        matrix that will be transformed into an orthonormal matrix Q and an upper triangular matrix R.
        
    Methods
    -------
    triangle_operation:
        Finds the householer matrix and uses it for the stepwise transformation of A into R and Q.
    triangle
        Returns the upper triangular matrix after Householder Transformation.
    vector 
        Returns the vectors of the Housholder matrix.
    """
    
    def __init__(self, A):
        """Matrix A"""
        
        self.A = A
        
    def triangle_operation(self):
        """Find the householer matrix and use it for the stepwise transformation of A into R.
        
        Create the vectors a, from the columns of matrix A, and e, from the columns of the 
        identity matrix I, the scalar sigma := ||a|| and the vector v for the 
        initialization of the Householder transformation using the class Reflection. 
            
        Returns
        -------
        R: numpy matrix
            The upper triangular matrix after Householder Transformation.
            
        list_v: list
            The vectors v for each step of the transformation."""
        
        shape = np.shape(self.A)
        m,n = shape
        
        R = np.copy(self.A)
        
        #Dependent on the shape of the matrix you have to do the transformation on a 
        #different number r of columns
        if m > n:
            r = n
        elif n > m:
            r = m
        else:
            r = n - 1
        
        # Create identity matrix I of same size as A
        I = np.zeros(m*r).reshape(m,r)
        
        I[:r] = np.identity(r)

        # Create list_v 
        list_v = []

        # write out vectors a and e of decreasing size from the columns of R and I 
        
        for j in list(range(r)): 
            a = [row[j] for row in R[j:]] # j'th column of A but only the m-i last rows.
            e = [row[j] for row in I[j:]] # same for the identity matrix
            
            a = np.array(a)
            e = np.array(e)
            sigma = np.linalg.norm(a) # this is the norm of the vector/column of A 
            v = a.reshape(m-j,1) + (np.dot(sigma, e.reshape(m-j,1))) # v = x + sigma * e
            list_v.append(v)

            H = Reflection(list_v[j]) # calculate the Housholder transformation for the vector v
            R = H * R # apply the transformation to the matrix A and obtain R stepwise

        
        return(R, list_v)
    
    def triangle(self):
        """Returns the upper triangular matrix after Householder Transformation.
        
        Returns
        -------
        R: numpy matrix
            the upper triangular matrix after Householder Transformation."""
        
        R = Householder.triangle_operation(self)[0]       
        
        return(R.round(10))
    
    def vector(self):
        """Returns the list of vectors v in the Housholder transformation.
        
        Retruns
        -------
         v_list: list
            The vectors v for each step of the transformation."""
        
        v_list = Householder.triangle_operation(self)[1]
        
        return(v_list)
    
    def obtain_Q(self):
        """Calculates the upper triangular matrix Q after QA decomposition.
        
        Use the list of housholder vectors to get the Householdermatrix H for 
        each step of
        the QR transformation,fill in the identitymatrix to match the
        shape of A, take te transpose and multiply to obtain Q.'
        
        Returns
        -------
        Q: numpy matrix
            The upper triangular matrix after QR decomposition using Householder 
            transformation.
        """
        
        #create the initial triangular matrix as a copy of the m x n - matrix A
            
        v_list = Householder.vector(self)
        n_n = len(v_list)
        q_m = len(v_list[0])
        
        H_list = []
        for i in list(range(n_n)):
            
            gamma = ((np.linalg.norm(v_list[i]))**2)/2
            vvtrans = v_list[i] * np.transpose(v_list[i])
            H =  np.identity((q_m-i)) - (vvtrans/gamma)
            
            print(H.shape)

            m_H, n_H = H.shape
            if m_H <  q_m:
                I = np.identity(q_m)
                x = y = i
                I [ x:x+H.shape[0], y:y+H.shape[1]] = H
                H = I
            H_list.append(H)
        
       # The transpose of Q is the result of the dot product H(n-1)...H1 
        
        len_H = len(H_list)

        H_temp = H_list[-1]
        
        for i in np.arange(len_H-1,0,-1):
            
            H_temp = np.matmul(H_temp, H_list[i-1])
        
        Q = np.transpose(H_temp)
            
        return(Q)


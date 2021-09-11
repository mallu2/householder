import sys
import os
import numpy as np

class Reflection:
    """
    Class for calculating the reflection of a, through the orthogonal matrix of v. 
    
    ...
    
    Attributes
    ----------
    v: numpy array 
        vector, shape(m,1), to initialize the Housholder transformation matrix.
    
    Methods
    -------
    __mul__(a)
        Transforms the vector a using the Housholder transformation 
        H = I - 2((vv(T)) / |v|**2).
    """
    
    def __init__(self, v): 
        """
        Vector v shape(m,1), to initialize the Housholder transformation.
        """
        self.v = v 
        
    def __mul__(self, a):
        """  
        Transforms the vector a using the Housholder transformation H = I - 2((vv(T)) / |v|**2).
        
        Calculate the matrix product for the 
        reflection of vector a orthogonal matrix H.
        The dimension of a needs to be bigger or equal 
        to the dimesion of v.
        If the dimension of v is smaller than the dimension of a, 
        v can be filled with zeros.
        
        Parameters
        ----------
        a: numpy array
            vector, shape(m,1) input vector to be transformed.
            
        Returns
        -------
        reflection: numpy array
            vector, shape(m,1), a reflected through H.
        """
        
        if len(a) < len(self.v):
            raise Exception('Length of a is smaller than length of v!')
        elif len(a) > len(self.v):
            differenceav = (len(a) - len(self.v))
            add = np.zeros(differenceav)
            velong = np.insert(self.v,0,add).reshape(len(a),1) # add zeros at the end
            self.v = velong
        else:
            self.v = self.v
            
        gamma = ((np.linalg.norm(self.v))**2)/2
        vvtrans = self.v * np.transpose(self.v)
        H =  np.identity((len(a))) - (vvtrans/gamma)
        reflection = np.dot(H,a)
        
        return(reflection) 

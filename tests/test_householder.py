#!/usr/bin/python3
# test_householder.py

import sys, os
import numpy as np

# import the packages that should be tested

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from package import Householder

class matrices:
    """This is a class that contains a number of matrices that can be used to test programs that to matrix calculations.
    
    Attributes
    ----------
    mn_33: numpy matrix
        3 x 3 matrix
    mn_22: numpy matrix
        2 x 2 matrix
    mn_44: numpy matrix
        4 x 4 matrix
    mn_43: numpy matrix
        4 x 3 matrix
    mn_34: numpy matrix
        3 x 4 matrix

    Methods
    -------
    m_33: numpy matrix
        returns 3 x 3 matrix
    m_22: numpy matrix
        returns 2 x 2 matrix
    m_44: numpy matrix
        returns 4 x 4 matrix
    m_43: numpy matrix
        returns 4 x 3 matrix
    m_34: numpy matrix
        returns 3 x 4 matrix"""

class matrices:
    """This is a class that contains a number of matrices that can be used to test programs that to matrix calculations.
    
    Attributes
    ----------
    mn_33: numpy matrix
        3 x 3 matrix
    mn_22: numpy matrix
        2 x 2 matrix
    mn_44: numpy matrix
        4 x 4 matrix
    mn_43: numpy matrix
        4 x 3 matrix
    mn_34: numpy matrix
        3 x 4 matrix

    Methods
    -------
    m_33: numpy matrix
        returns 3 x 3 matrix
    m_22: numpy matrix
        returns 2 x 2 matrix
    m_44: numpy matrix
        returns 4 x 4 matrix
    m_43: numpy matrix
        returns 4 x 3 matrix
    m_34: numpy matrix
        returns 3 x 4 matrix"""
    
    
    
    def __init__(self):
        self.mn_33 = np.array([[1,2,5],[3,4,7],[2,9,7]])
        self.mn_22 = np.array([[1,2],[3,4]])
        self.mn_44 = np.array([[1,2,5,7],[3,4,5,5],[4,1,4,8],[5,8,9,1]])

        self.mn_43 = np.array([2,3,4,(-2),4,7,5,9,1,3,5,6])
        self.mn_43 = self.mn_43.reshape(4,3)
        self.mn_34 = self.mn_43.reshape(3,4)

    def  m_33(self):
        "Returns 3 x 3 matrix."
        return(self.mn_33)

    def  m_22(self):
        "Returns 2 x 2 matrix."
        return(self.mn_22)

    def  m_44(self):
        "Returns 4 x 4 matrix."
        return(self.mn_44)

    def  m_43(self):
        "Returns 4 x 3 matrix."
        return(self.mn_43)

    def  m_34(self):
        "Returns 3 x 4 matrix."
        return(self.mn_34)


def QR_decomb_t(A):
    """This function tests QR decompositions of a matrix A by the Householder class.
    
    Parameters
    ----------
    A : numpy matrix
        input matrix that is decomposed.
    
    Returns
    -------
    result == A: array of bool
        should be all True of the algorith is working correctly."""
    
    At = A.copy()
    R = Householder.triangle(At)
    v = Householder.vector(At)
    Q = Householder.obtain_Q(At)
    result = np.dot(Q, R)
    result = np.rint(result).astype(int)
    
    print ("If all fields are true you are good!",result == A)  


def test():
    "Check if the QR desomposition of A results in QR = A."
    for matrix in [matrices().m_22(), matrices().m_33(), matrices().m_34(), matrices().m_43(), matrices().m_44()]:
        QR_decomb_t(matrix)

test()
     

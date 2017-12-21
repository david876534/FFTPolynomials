#   This project's goal is to calculate the convolution of two vector represented polynomials - 
#   and to compare the speed of the Diverse Fourier Transform and the Fast Fourier Transform.
#   The speed up of the nlog(n) algorithm over the quadratic time DFT is evident at n = 16
#-----------------------------------Running Details-------------------------------------
#   Run this program with the command python3 fftPolynomialMult.py
#   Input can be entered at the prompt, the program terminates after a single run.

import math
import random
import sys
import time

sys.setrecursionlimit(10000000)

from random import *

NGlobal = 0

def DFT(A):
    n = len(A)
    return [sum(A[k] * calcOmega(-i * k) for k in range(n))
            for i in range(n)]

def DFTINV(A):
    n = len(A)
    return [round((sum(A[k] * calcOmega(i * k) for k in range(n)) / n).real)
            for i in range(n)]

def FFT(A):
    global NGlobal
    NGlobal = len(A)
    omega = math.cos(2*math.pi/NGlobal) + math.sin(2*math.pi/NGlobal)*1j
    return fft(A, omega, 1) 

def FFTINV(A):
    v = FFT(A)
    ans = []
    ans.append(round((v[0] / NGlobal).real))
    for i in range(len(v) - 1, 0, -1):
        ans.append(round((v[i] / NGlobal).real))
    return ans

def fft(A, omega, n):
    if calcOmega(n).real == 1:
        return A
    s = fft(evenHalf(A), calcOmega(n * 2), n * 2)
    sp = fft(oddHalf(A), calcOmega(n * 2), n * 2)
    R = []
    for j in range(0, len(A)):
        R.append(0 + 0j)
    for i in range(0, int(len(A) / 2)):
        R[i] = s[i] + sp[i] * omega**i
        R[i + int(len(A)/2)] = s[i] - sp[i] * omega**i
    return R

def calcOmega(n) :
    a = (math.e**(2 * math.pi * 1j * n/NGlobal))
    return a

def evenHalf(A):
    even = []
    for i in range(0, int(len(A) / 2)):
        even.append(A[i * 2])
    return even

def oddHalf(A):
    odd = []
    for i in range(0, int(len(A) / 2)):
        odd.append(A[i*2 + 1])
    return odd

def fillRand(n):
    A = []
    for a in range(0, n):
        A.append(randint(1, 1000) % 2 + 1)
    return A

def pad(A):
    for i in range(0, len(A)):
        A.append(0)
    return A

def compMultiply(A, B):
    R = []
    for i in range(0, len(A)):
        R.append(A[i] * B[i])
    return R

def getReal(A):
    RET = []
    for i in range(len(A)):
        RET.append(A[i].real)
    return RET

def main():
    n = int(input("Please enter an n value (must be a power of 2): "))
    A = []
    for i in range(n):
        A.append(int(input("enter a {} degree vector item for A: ".format(i))))
    B = []
    for i in range(n):
        B.append(int(input("enter a {} degree vector item for B: ".format(i))))
      
    print("Vector 1: ")
    print(A)
    print("Vector 2: ")
    print(B)
    A = pad(A)
    B = pad(B)

    startFFT = time.time()
 
    R = FFT(A)
    S = FFT(B)


    ans = FFTINV(compMultiply(R, S))
    
    endFFT = time.time()
    print("Time to run FFT: {}".format(endFFT-startFFT))
    print("Convolution: {}".format(ans))
        
    
    startDFT = time.time()
    R = DFT(A)
    S = DFT(B)
    
    ans = DFTINV(compMultiply(R, S))
    endDFT = time.time()
    print("Time to run DFT: {}".format(endDFT-startDFT))
    print("Convolution: {}".format(ans))
    

if __name__ == '__main__':
    main()

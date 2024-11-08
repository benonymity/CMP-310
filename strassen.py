# 
# By Benjamin Bassett (benonymity) on 11.7.24
# 
# For CMP-310 with Dr. Serang at Hillsdale College
# https://github.com/benonymity/CMP-310
# 

import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt

# Implemented from scratch to make sure the @ operator is doing what I think it is
def naive(A: np.ndarray, B: np.ndarray):
    n = A.shape[0]
    result = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                result[i,j] += A[i,k] * B[k,j]
    return result

def strassen(A: np.ndarray, B: np.ndarray, leaf_size: int = 64):
    # How it works:
    #
    # The result matrix C is computed as:
    #
    #     C     =     Result
    # |C11 C12| |P5+P6-P1+P4    P1+P2    |
    # |C21 C22| |   P3+P4    P5+P7+P2-P3 |
    #
    # Where P1 through P7 are:
    # P1 = (A11 + A12)B22
    # P2 = (B12 - B22)A11
    # P3 = (A21 + A22)B11
    # P4 = (B21 - B11)A22
    # P5 = (A11 + A22)(B11 + B22)
    # P6 = (A12 - A22)(B21 + B22)
    # P7 = (A21 - A11)(B11 + B12)

    # Won't check for squareness or do padding because we're spoon-feeding nice inputs here
    n = A.shape[0]

    if n <= leaf_size:
        return A @ B # Banged my head against this for a while——turns out A * B ≠ A @ B...

    # Split the matrices into quadrants
    mid = n//2
    A11 = A[:mid, :mid]
    A12 = A[:mid, mid:] # This matrix notation seems backwards——why are the rows indexed first??
    A21 = A[mid:, :mid]
    A22 = A[mid:, mid:]

    B11 = B[:mid, :mid]
    B12 = B[:mid, mid:]
    B21 = B[mid:, :mid]
    B22 = B[mid:, mid:]

    # Compute the P's
    P1 = strassen((A11 + A12), B22, leaf_size)
    P2 = strassen(A11, (B12 - B22), leaf_size) # Why isn't this commutative?
    P3 = strassen((A21 + A22), B11, leaf_size)
    P4 = strassen(A22, (B21 - B11), leaf_size) # Same here
    P5 = strassen((A11 + A22), (B11 + B22), leaf_size) # The magic happens!
    P6 = strassen((A12 - A22), (B21 + B22), leaf_size)
    P7 = strassen((A21 - A11), (B11 + B12), leaf_size)

    C11 = P5 + P6 + P4 - P1
    C12 = P1 + P2
    C21 = P3 + P4
    C22 = P5 + P7 + P2 - P3

    C = np.vstack((np.hstack((C11, C12)), np.hstack((C21, C22)))) # Combine the quadrants with numpy magic
    return C

if __name__ == "__main__":
    sys.setrecursionlimit(10000) # Just in case ;)
    
    sizes = [32, 64, 128, 256, 512, 1024] # Try a bunch of powers of two to avoid padding
    naive_times = []
    strassen_times = []
    leaf_size = 128 # This can be optimized somehow

    for n in sizes:
        print(f"\nTesting matrices of size {n}x{n}")
        A = np.random.randint(0, 10, (n, n), dtype=np.int64)
        B = np.random.randint(0, 10, (n, n), dtype=np.int64)
        
        # Time naive
        t1 = time()
        naive_result = A @ B
        t2 = time()
        naive_time = t2 - t1
        print(f"Naive multiplication took {naive_time:.2f} seconds")
        naive_times.append(naive_time)
        
        # Time Strassen
        t1 = time()
        strassen_result = strassen(A, B)
        t2 = time()
        strassen_time = t2 - t1
        print(f"Strassen took {strassen_time:.2f} seconds")
        strassen_times.append(strassen_time)
        
        # Verify the results match
        if not np.allclose(naive_result, strassen_result):
            print("Results don't match :(")
            print(f"Naive result: {naive_result}")
            print(f"Strassen result: {strassen_result}")
            break

    # Make a log-log plot of the runtimes!
    plt.plot(sizes, naive_times, label="Naive")
    plt.plot(sizes, strassen_times, label="Strassen")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig(f"strassen_comparison_{leaf_size}_leaves.png")

    print(f"Saved plot to strassen_comparison_{leaf_size}_leaves.png")
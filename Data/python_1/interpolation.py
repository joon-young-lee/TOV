import numpy as np

def my_cubic_spline_interpolation(x, fx, xi): # len(x) - 1 splines
    # S_i(x) = fx + b_i * (x-xi) + c_i * (x-xi) ** 2 + d_i * (x-xi) ** 3
    # Each spline has 3 unknown coefficients (we alrady know the constant a_i)
    if len(x) != len(fx):
        return 'error' 
    elif max(xi) > max(x) or min(xi) < min(x):
        return 'error'
    n = len(x)
    N = 3 * (n-1) 
    A = np.zeros((N, N))
    b = np.zeros(N)
    a = fx # constant a
    # index for b_i 0, 3, 6, 9 ...
    # index for c_i 1, 4, 7, 10 ...
    # index for d_i 2, 5, 8, 11 ...
    for k in range(N):
        if 0 <= k <= n-2:
                b[k] = fx[k+1] - fx[k]
    for i in range(N):
        # 1st rows for equal fx
        if i in range(n-1): # 0 <= i <= n-2 (#: n-1)
            for l in range(3):
                A[i, 3 * i + l] = (x[i+1] - x[i]) ** (l + 1)
        
        # 2nd rows for equal deriv
        elif i in range(n-1, 2 * n - 3): # n-1 <= i <= 2 * n - 5 (#: n-2)
            for l in range(4):
                if l == 0:
                    A[i, 3 * (i - (n-1)) + l] = 1
                elif l == 1:
                    A[i, 3 * (i - (n-1)) + l] = 2 * (x[i - (n-2) + 1] - x[i - (n-2)])
                elif l == 2:
                    A[i, 3 * (i - (n-1)) + l] = 3 * (x[i - (n-2) + 1] - x[i - (n-2)]) ** 2
                elif l == 3:
                    A[i, 3 * (i - (n-1)) + l] = -1
        
        # 3rd row for equal second deriv
        elif i in range(2 * n - 3, 3 * n - 5): # 2 * n - 4 <= i <= 3 * n - 7 (#: n-2)
            for l in range(5):
                if l == 1:
                    A[i, 3 * (i- (2*n-3)) + l] = 1
                elif l == 4:
                    A[i, 3 * (i- (2*n-3)) + l] = -1
                elif l == 2:
                    A[i, 3 * (i- (2*n-3)) + l] = 3 * (x[i - 2*(n-2) + 1] - x[i - 2*(n-2)])
        
        # Boundaries
        elif i == 3 * n - 5: #(#: 1)
            A[i, 1] = 1
        
        elif i == 3 * n - 4: #(#: 1)
            A[i, N - 2] = 2
            A[i, N - 1] = 6
    
    coe = np.linalg.solve(A, b)
    # we have coefficients
    # now we have to get the right coefficients for xi
    H = len(xi) # number of points we want to interpolate
    fxi = np.zeros(H)
    for i in range(H):
        for j in range(n - 1):
            if xi[i] == x[j]:
                fxi[i] = fx[j]
            elif xi[i] == x[j+1]:
                fxi[i] = fx[j+1]
            elif x[j] < xi[i] < x[j+1]:
                fxi[i] = a[j] + \
                    coe[3 * j] * (xi[i] - x[j]) + \
                    coe[3 * j + 1] * (xi[i] - x[j]) ** 2 + \
                    coe[3 * j + 2] * (xi[i] - x[j]) ** 3
    
    return fxi

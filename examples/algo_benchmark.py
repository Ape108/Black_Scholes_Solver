import time
import random

def gaussian_elimination(A, d):
    """
    Naive O(N^3) Gaussian Elimination for a dense matrix.
    """
    n = len(A)
    
    # Forward elimination
    for i in range(n):
        pivot = A[i][i]
        if pivot == 0.0:
            raise ValueError("Zero pivot encountered.")
            
        for j in range(i + 1, n):
            if A[j][i] != 0.0:
                factor = A[j][i] / pivot
                for k in range(i, n):
                    A[j][k] -= factor * A[i][k]
                d[j] -= factor * d[i]
                
    # Back substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        row_sum = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (d[i] - row_sum) / A[i][i]
        
    return x

def thomas_algorithm(a, b, c, d):
    """
    O(N) Thomas Algorithm for a tridiagonal matrix.
    a: main diagonal
    b: upper diagonal
    c: lower diagonal
    """
    n = len(d)
    
    # Create copies to avoid mutating the original arrays
    a_mod = a[:]
    d_mod = d[:]
    
    # 1. Forward Sweep
    for i in range(1, n):
        factor = c[i-1] / a_mod[i-1]
        a_mod[i] -= factor * b[i-1]
        d_mod[i] -= factor * d_mod[i-1]
        
    # 2. Backward Substitution
    x = [0.0] * n
    x[-1] = d_mod[-1] / a_mod[-1]
    
    for i in range(n - 2, -1, -1):
        x[i] = (d_mod[i] - b[i] * x[i+1]) / a_mod[i]
        
    return x

def run_benchmarks():
    # Test different grid sizes (N)
    sizes = [100, 250, 500, 1000, 2000, 5000]
    
    print(f"{'Matrix Size (N)':<15} | {'O(N^3) Gaussian':<15} | {'O(N) Thomas':<15} | {'Speedup':<10}")
    print("-" * 65)
    
    for n in sizes:
        # 1. Generate random tridiagonal data
        a_diag = [random.uniform(5.0, 10.0) for _ in range(n)] # Main
        b_diag = [random.uniform(1.0, 3.0) for _ in range(n-1)] # Upper
        c_diag = [random.uniform(1.0, 3.0) for _ in range(n-1)] # Lower
        d_vec = [random.uniform(10.0, 50.0) for _ in range(n)]  # RHS
        
        # 2. Reconstruct the dense matrix A for Gaussian Elimination
        A_dense = [[0.0] * n for _ in range(n)]
        for i in range(n):
            A_dense[i][i] = a_diag[i]
            if i < n - 1:
                A_dense[i][i+1] = b_diag[i]
            if i > 0:
                A_dense[i][i-1] = c_diag[i-1]
                
        # 3. Time Gaussian Elimination
        start_gauss = time.perf_counter()
        # Pass copies so we don't destroy the original lists
        x_gauss = gaussian_elimination([row[:] for row in A_dense], d_vec[:])
        gauss_time = time.perf_counter() - start_gauss
        
        # 4. Time Thomas Algorithm
        start_thomas = time.perf_counter()
        x_thomas = thomas_algorithm(a_diag, b_diag, c_diag, d_vec)
        thomas_time = time.perf_counter() - start_thomas
        
        # 5. Output results
        speedup = gauss_time / thomas_time if thomas_time > 0 else float('inf')
        
        print(f"{n:<15} | {gauss_time:<15.5f} | {thomas_time:<15.5f} | {speedup:.1f}x")

if __name__ == "__main__":
    run_benchmarks()
# Numerical Programming Library
**Author**: 22400807  
**Course**: Numerical Programming — Handong Global University  
**Language**: C++ (MSVS 2022)  
**Files**: `myNP_22400807.h/cpp`, `myMatrix_22400807_2026.h/cpp`

---

## Project Structure

```
include/
├── myMatrix_22400807_2026.h   # Matrix struct & operations
└── myNP_22400807.h            # Numerical methods

src/
├── myMatrix_22400807_2026.cpp
└── myNP_22400807.cpp
```

---

## myMatrix — Matrix Operations

### Data Structure
```cpp
typedef struct {
    double** at;   // 2D array: access via A.at[i][j]
    int rows, cols;
} Matrix;
```

### Basic Operations
| Function | Description |
|---|---|
| `zeros(rows, cols)` | Create zero matrix |
| `ones(rows, cols)` | Create ones matrix |
| `eye(rows, cols)` | Create identity matrix |
| `createMat(rows, cols)` | Allocate matrix |
| `freeMat(A)` | Free matrix memory |
| `initMat(A, val)` | Fill matrix with value |
| `copyMatrix(Aout, Bin)` | Copy B into A |
| `arr2Mat(array, rows, cols)` | 1D array → Matrix |
| `txt2Mat(path, name)` | Text file → Matrix |
| `printMat(A, "name")` | Print matrix |

### Matrix Arithmetic
| Function | Description |
|---|---|
| `addMat(A, B)` | A + B |
| `subMat(A, B)` | A - B |
| `multMat(A, B)` | A × B |
| `smultMat(A, k)` | k × A (scalar multiply) |
| `transpose(A)` | Aᵀ |
| `invMat(A)` | A⁻¹ |
| `norm(v)` | Vector 2-norm |

### Linear System Solvers
| Function | Description |
|---|---|
| `solveLinear(A, b)` | Solve Ax = b |
| `fwdSub(L, b)` | Forward substitution |
| `bwdSub(U, y)` | Backward substitution |
| `GaussElimination(A, b, L, U, B)` | Gauss elimination |
| `LUdecomp(L, U, A)` | LU decomposition |
| `LUdecomp_pivot(L, U, P)` | LU with pivoting |

### Eigenvalue & Eigenvector
| Function | Description |
|---|---|
| `eigval_student(A)` | Eigenvalues (QR iteration) |
| `eigvec_student(A)` | Eigenvectors |
| `eig(V, D, A)` | Both eigenvalues & eigenvectors |
| `QRdecomp_student(Q, R, A)` | QR decomposition |

```cpp
// Usage
Matrix A = arr2Mat(data, 3, 3);
Matrix V = zeros(3, 3);   // eigenvectors
Matrix D = zeros(3, 3);   // eigenvalues (diagonal)
eig(V, D, A);
// D.at[i][i] = i-th eigenvalue (sorted descending)
// V.at[:,i]  = i-th eigenvector
```

### Curve Fitting (Matrix version)
| Function | Description |
|---|---|
| `linearFit_mat(X, Y)` | Linear regression → [a0, a1] |
| `polyFit_mat(X, Y, n)` | Polynomial regression (nth order) |
| `expFit_mat(X, Y)` | Exponential regression y=C0*e^(C1*x) |
| `multiFit_mat(X, Y)` | **Multivariate** linear regression |

```cpp
// multiFit_mat — Multivariate Linear Regression
// Model: Y = a0 + a1*X1 + a2*X2 + ... + aD*XD
// W* = (A^T A)^-1 A^T Y
// X: m x D,  Y: m x 1,  W: (D+1) x 1

Matrix X = arr2Mat(xdata, m, D);   // D variables
Matrix Y = arr2Mat(ydata, m, 1);
Matrix W = multiFit_mat(X, Y);
// W.at[0][0] = a0, W.at[1][0] = a1, ...
```

### Normalization (for Gradient Descent)
| Function | Description |
|---|---|
| `normalizeMat(A)` | Normalize to [0, 1] |
| `unnormMat(Z, X_raw)` | Restore original scale |
| `maxMat(A)` | Max value in matrix |
| `minMat(A)` | Min value in matrix |

---

## myNP — Numerical Methods

### Nonlinear Equation Solvers
| Function | Signature |
|---|---|
| `bisection` | `bisection(func, a, b, tol)` |
| `newtonRaphson` | `newtonRaphson(func, dfunc, x0, tol)` |

### Interpolation
| Function | Signature |
|---|---|
| `linearInterp` | `linearInterp(x_arr, y_arr, n, xq)` |

### Differentiation
| Function | Description |
|---|---|
| `gradient1D(dydx, x, y, m)` | 1st derivative from data |
| `gradientFunc(dydx, func, x, m)` | 1st derivative from function |
| `acceleration(dy2dx2, x, y, m)` | 2nd derivative from data |

### Integration
| Function | Description |
|---|---|
| `simpson(x, y, m)` | Simpson's rule from data |
| `integral(func, a, b, n)` | Simpson's rule from function |

### ODE Solvers

#### 1st Order ODE: `dy/dt = f(t, y)`
| Function | Method |
|---|---|
| `odeEU(y, func, t0, tf, h, y0)` | Euler Explicit |
| `odeRK2(y, func, t0, tf, h, y0)` | Runge-Kutta 2nd (Modified Euler) |
| `odeRK3(y, func, t0, tf, h, y0)` | Runge-Kutta 3rd |

```cpp
// 1st order ODE example: dy/dt = -2y
double myODE(const double t, const double y) {
    return -2 * y;
}
int N = (tf - t0) / h + 1;
double y[N] = {0};
odeRK2(y, myODE, 0, 5, 0.1, 1.0);
```

#### 2nd Order ODE: `y'' + cy' + ky = u(t)`
Convert to system: `dy/dt = v`, `dv/dt = -ky - cv + u`

| Function | Method |
|---|---|
| `sys2EU(y, v, func, t0, tf, h, y0, v0)` | Euler Explicit |
| `sys2RK2(y, v, func, t0, tf, h, y0, v0)` | Runge-Kutta 2nd |
| `sys2RK3(y, v, func, t0, tf, h, y0, v0)` | Runge-Kutta 3rd |

```cpp
// 2nd order ODE: y'' + 0.5y' + y = 0
// k=1, c=0.5, u=0
void myODESys(double dYdt[], const double t, const double Y[]) {
    dYdt[0] = Y[1];                          // dy/dt = v
    dYdt[1] = -1*Y[0] - 0.5*Y[1] + 0;      // dv/dt = -ky - cv + u
}
int N = (tf - t0) / h + 1;
double y[N] = {0}, v[N] = {0};
sys2RK2(y, v, myODESys, 0, 10, 0.01, 1.0, 0.0);
//             y0=1, v0=0 (initial displacement, velocity)
```

### Linear Regression (Array version)
| Function | Description |
|---|---|
| `linearRegression(z, x, y, n)` | Linear: y = a0 + a1*x |
| `polyFit(z, x, y, n, order)` | Polynomial nth order |
| `expRegression(z, x, y, n)` | Exponential: y = C0*e^(C1*x) |

### Gradient Descent

```cpp
// Data Fitting via Gradient Descent
// Model: Y = a0 + a1*X1 + ... + aD*XD

// Step 1. Normalize
Matrix X = normalizeMat(X_raw);   // X: M x D
Matrix W = zeros(D+1, 1);         // W size rule: zeros(D+1, 1)
                                   // 1D -> zeros(2,1)
                                   // 2D -> zeros(3,1)
// Step 2. Run
gradDescFit(lossFunc_linear, lossGradFunc_linear, W, X, Y, eta, itrMax);

// Step 3. Unnormalize
Matrix W_est = unnormMat(W, X_raw);
```

| Function | Description |
|---|---|
| `lossFunc_linear(W, X, Y)` | MSE Loss: (1/M)Σ(y-ŷ)² |
| `lossGradFunc_linear(W, X, Y)` | Gradient: -2/M * Σ(r·A) |
| `gradDescFit(lossFunc, lossGradFunc, W, X, Y, eta, itrMax)` | GD solver |

---

## Quick Reference — W Size Rule

```
Model                  X shape    W size
─────────────────────────────────────────
y = a0 + a1*x          m x 1     zeros(2,1)
y = a0 + a1*x1 + a2*x2 m x 2    zeros(3,1)
y = a0 + ... + aD*xD   m x D     zeros(D+1,1)
```

---

## Memory Management

**반드시** 사용 후 `freeMat()` 호출:
```cpp
Matrix A = zeros(3, 3);
// ... use A ...
freeMat(A);   // must free!
```

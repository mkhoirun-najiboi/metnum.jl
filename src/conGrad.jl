"""
    conGrad(A,B,Xawal)
adalah fungsi yang digunakan untuk mencari solusi SPL `Ax=B` secara iteratif dengan
nilai tebakan awal `Xawal` menggunakan metode conjugate gradient.

# Example
```jl
julia> A = [2  3  3
        3  5  4
        3  4  6];

julia> B = [3,4,6];

julia> X,flag,err,M = conGrad(A,B,[0 0 0.9]);

julia> X
3-element Array{Float64,1}:
  2.380760976839332e-6
 -1.220333023839802e-6
  0.9999996142293243

julia> flag
0

julia> err
8.019009497939746e-8

julia> M
229×5 Array{Float64,2}:
   0.0  0.0          0.0          0.9       NaN
   1.0  0.0266764    0.0355685    0.953353    0.0694497
   2.0  0.0154769    0.000266424  0.975767    0.0432903
   3.0  0.0193622    0.00858077   0.976843    0.00924021
   ⋮
 225.0  5.10651e-6  -1.83459e-6   0.999999    1.32762e-6
 226.0  2.29335e-6  -1.29137e-6   1.0         3.08384e-6
 227.0  2.4522e-6   -1.23209e-6   1.0         1.73709e-7
 228.0  2.38076e-6  -1.22033e-6   1.0         8.01901e-8
```
return solusi `X` dengan `flag` bernilai 0 jika metode conjugate gradient
berhasil menemukan solusi dan gagal jika tidak 0.
Serta, matriks `M` yang berisi catatan proses tiap iterasi `[k, X', error]`
"""
function conGrad(A,B,Xawal)
    if ~(isposdef(A) && issymmetric(A))
        error("matriks A harus simetrik definit positif")
    end
    delta = 10^-7;
    maxi = 1000;
    flag = 1;
    X = Xawal[:];
    r = B-A*X
    d = r
    M = [0 X' NaN];
    for k = 1:maxi
        Xlama = X
        rlama = r
        dlama = d
        X = Xlama + ((rlama'rlama)/(dlama'A*dlama))*dlama
        r = rlama - ((rlama'rlama)/(dlama'A*dlama))*A*dlama
        d = r - ((r'r)/(rlama'rlama))*dlama
        err = norm(X-Xlama)
        M = [M; [k X' err] ];
        if err<delta
            flag = 0;
            break
        end
    end
    err = M[end,end]
    return X, flag, err, M
end

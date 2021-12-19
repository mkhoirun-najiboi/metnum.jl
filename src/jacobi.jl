"""
    jacobi(A, B, Xawal)
adalah fungsi yang digunakan untuk mencari solusi SPL `Ax=B` secara iteratif dengan
nilai tebakan awal `Xawal` menggunakan metode jacobi.

# Example
```jl
julia> A = [4 -1 1
            4 -8 1
           -2 1 5];

julia> B = [7;-21;15];

julia> Xa = [1,2,2];

julia> X,flag,err,M = jacobi(A,B,Xa);

julia> X
3-element Array{Float64,1}:
 1.9999999884128572
 3.9999999907302857
 3.0000000018539428

julia> flag
0

julia> err
3.777053140308286e-8

julia> M
18×5 Array{Float64,2}:
  0.0  1.0      2.0    2.0     NaN
  1.0  1.75     3.375  3.0       1.85826
  2.0  1.84375  3.875  3.025     0.509327
  3.0  1.9625   3.925  2.9625    0.143205
  ⋮
 14.0  2.0      4.0    3.0       1.00721e-6
 15.0  2.0      4.0    3.0       2.83194e-7
 16.0  2.0      4.0    3.0       1.37804e-7
 17.0  2.0      4.0    3.0       3.77705e-8
```
return solusi `X` dengan `flag` bernilai 0 jika metode jacobi berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, X', error]`
"""
function jacobi(A, B, Xawal)
    delta = 10^-7;
    maxi = 100;
    flag = 1;
    D = Diagonal(diag(A))
    R = A - D;
    X = Xawal[:];
    M = [0 X' NaN];
    for k = 1:maxi
        Xlama = X;
        X = inv(D)*(B - R*Xlama);
        err = norm(X-Xlama);
        M = [M; [k X' err] ];
        if err<delta || norm(B-A*X)<delta
            flag = 0;
            break
        end
    end
    err = M[end,end]
    return X, flag, err, M
end

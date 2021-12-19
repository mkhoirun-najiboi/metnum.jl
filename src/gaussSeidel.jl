"""
    gaussSeidel(A,B,Xawal)
adalah fungsi yang digunakan untuk mencari solusi SPL `Ax=B` secara iteratif dengan
nilai tebakan awal `Xawal` menggunakan metode Gauss-Seidel.

# Example
```jl
julia> A = [4 -1 1
            4 -8 1
           -2 1 5];

julia> B = [7;-21;15];

julia> Xa = [1,2,2];

julia> X,flag,err,M = gaussSeidel(A,B,Xa);

julia> X
3-element Array{Float64,1}:
 1.99999999743314
 3.999999998137355
 2.999999999345785

julia> flag
0

julia> err
2.1378379298357413e-8

julia> M
11×5 Array{Float64,2}:
  0.0  1.0      2.0      2.0      NaN
  1.0  1.75     3.75     2.95       2.12779
  2.0  1.95     3.96875  2.98625    0.298606
  3.0  1.99562  3.99609  2.99903    0.0547054
  ⋮
  7.0  2.0      4.0      3.0        1.21519e-5
  8.0  2.0      4.0      3.0        1.34052e-6
  9.0  2.0      4.0      3.0        1.85294e-7
 10.0  2.0      4.0      3.0        2.13784e-8
```
return solusi `X` dengan `flag` bernilai 0 jika metode Gauss-Seidel berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, X', error]`
"""
function gaussSeidel(A,B,Xawal)
    delta = 10^-7;
    maxi = 100;
    flag = 1;
    D = Diagonal(diag(A))
    R = tril(A);
    U = triu(A,1);
    X = Xawal[:];
    M = [0 X' NaN];
    for k = 1:maxi
        Xlama = X;
        X = inv(R)*(B - U*Xlama);
        err = norm(X-Xlama);
        M = [M; [k X' err] ];
        if err<delta
            flag = 0;
            break
        end
    end
    err = M[end,end]
    return X, flag, err, M
end

"""
    rekons(A,B,Xawal)
adalah fungsi yang digunakan untuk mencari solusi SPL `Ax=B` secara iteratif dengan
nilai tebakan awal `Xawal` menggunakan metode rekonstruksi aljabar.

# Example
```jl
julia> A = [4 -1 1
            4 -8 1
           -2 1 5];

julia> B = [7;-21;15];

julia> Xa = [1,2,2];

julia> X,flag,err,M = rekons(A,B,Xa);

julia> X
3-element Array{Float64,1}:
 1.9999987931831584
 3.9999993707528194
 2.9999996431226994

julia> flag
0

julia> err
7.477270454885512e-8

julia> M
18×5 Array{Float64,2}:
  0.0  1.0       2.0      2.0      NaN
  1.0  0.677229  3.45151  2.58059    0.658544
  2.0  1.51207   3.75397  2.85403    0.0644652
  3.0  1.79762   3.89498  2.94005    0.0145717
  ⋮
 14.0  1.99998   3.99999  3.0        9.83924e-7
 15.0  1.99999   4.0      3.0        4.16762e-7
 16.0  2.0       4.0      3.0        1.76529e-7
 17.0  2.0       4.0      3.0        7.47727e-8
```
return solusi `X` dengan `flag` bernilai 0 jika metode rekonstruksi aljabar
berhasil menemukan solusi dan gagal jika tidak 0.
Serta, matriks `M` yang berisi catatan proses tiap iterasi `[k, X', error]`
"""
function rekons(A,B,Xawal)
    delta = 10^-7;
    maxi = 100;
    flag = 1;
    X = Xawal[:];
    Xlama = X;
    M = [0 X' NaN];
    for k = 1:maxi
        for i = 1:length(B)
            Xlama = X;
            X = Xlama + A[i,:]*(B[i]-A[i,:]'Xlama)/(A[i,:]'A[i,:])
        end
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

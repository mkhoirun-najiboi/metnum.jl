"""
    LUtanpaP(A)
adalah fungsi yang digunakan untuk memfaktorisasi matriks `A` sembarang menjadi
dua matriks, yaitu matriks segitiga atas `U` dan bawah `L` sedemikian sehingga `LU=A`

# Examples
```jl
julia> A = [ 2 5 4 2 3
            4 1 0.4 6 5
            0.9 4 8 7 10
            5 7 0.5 2 7
            2 6 8 9 4];

julia> L,U = LUtanpaP(A);

julia> L
5×5 Array{Float64,2}:
 1.0    0.0        0.0       0.0      0.0
 2.0    1.0        0.0       0.0      0.0
 0.45  -0.194444   1.0       0.0      0.0
 2.5    0.611111  -1.02824   1.0      0.0
 1.0   -0.111111   0.668235  1.17806  1.0

julia> U
5×5 Array{Float64,2}:
 2.0   5.0   4.0      2.0        3.0
 0.0  -9.0  -7.6      2.0       -1.0
 0.0   0.0   4.72222  6.48889    8.45556
 0.0   0.0   0.0      2.44988    8.80541
 0.0   0.0   0.0      0.0      -15.1347
```
return matriks `L` dan `U`.
"""
function LUtanpaP(A);
    # Definisikan matriks L sebagai pencatat pengali.
    n,n = size(A);
    L = zeros(n,n);
    # Lakukan operasi baris dasar terhadap matriks
    Aug = copy(A);
    for p = 1:n-1
        # Jika pivot bernilai nol, maka gagal.
        if Aug[p,p]==0
            error("Pivot bernilai nol");
        end
        # Lakukan eliminasi lalu simpan pengali pada matriks L.
        for i = p+1:n
            k = Aug[i,p]/Aug[p,p];
            Aug[i,1:n] = Aug[i,1:n] - k*Aug[p,1:n];
            L[i,p] = k;
        end
    end
    U = Aug;
    L = L .+ I(n);
    return L, U
end

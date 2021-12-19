"""
    LUdenganP(A)
adalah fungsi yang digunakan untuk memfaktorisasi matriks `A` sembarang menjadi
dua matriks, yaitu matriks segitiga atas `U` dan bawah `L` serta matriks permutasi P
sedemikian sehingga `LU=PA`.

# Examples
```jl
julia> A = [ 2 5 4 2 3
            4 1 0.4 6 5
            0.9 4 8 7 10
            5 7 0.5 2 7
            2 6 8 9 4];

julia> L,U,P = LUdenganP(A);

julia> L
5×5 Array{Float64,2}:
 1.0    0.0       0.0        0.0       0.0
 0.8    1.0       0.0        0.0       0.0
 0.18  -0.595652  1.0        0.0       0.0
 0.4   -0.695652  0.986094   1.0       0.0
 0.4   -0.478261  0.480405  -0.537685  1.0

julia> U
5×5 Array{Float64,2}:
  5.0           7.0  0.5   2.0       7.0
  0.0          -4.6  0.0   4.4      -0.6
  1.11022e-16   0.0  7.91  9.26087   8.38261
 -1.09478e-16   0.0  0.0   2.12879  -7.48343
 -1.122e-16     0.0  0.0   0.0      -8.13773

julia> P
5×5 Array{Bool,2}:
 0  0  0  1  0
 0  1  0  0  0
 0  0  1  0  0
 0  0  0  0  1
 1  0  0  0  0
```
return matriks `L`, `U`, dan `P`.
"""
function LUdenganP(A)
    # Definisikan matriks L dan P.
    n,n = size(A);
    L = zeros(n,n);
    P = Array(I(n));
    # Lakukan operasi baris dasar terhadap matriks gandeng
    Aug = copy(A)
    for p = 1:n-1
        # Lakukan operasi pindah baris untuk menentukan nilai pivot
        val,j = findmax(abs.(Aug[p:n,p]));
        # Pivoting U
        C = Aug[p,:];
        Aug[p,:] = Aug[j+p-1,:];
        Aug[j+p-1,:] = C;
        # Pivoting L
        C = L[p,:];
        L[p,:] = L[j+p-1,:];
        L[j+p-1,:] = C;
        # Pivoting P
        C = P[p,:];
        P[p,:] = P[j+p-1,:];
        P[j+p-1,:] = C;
        # Jika pivot bernilai nol, maka gagal.
        if Aug[p,p]==0
            error("Matriks Singular")
        end
        # Lakukan eliminasi dan simpan pengali.
        for i = p+1:n
            k = Aug[i,p]/Aug[p,p]
            Aug[i,1:n] = Aug[i,1:n] - k*Aug[p,1:n];
            L[i,p] = k;
        end
    end
    U = Aug;
    L = L + I(n);
    return L,U,P
end

"""
    elimGaussNonPivoting(A,B)
adalah fungsi untuk menyelesaikan SPL `Ax=B` dengan `A` adalah matriks sembarang
dan `B` adalah vektor ruas kanan menggunakan metode eliminasi Gauss tanpa pivoting.

# Examples
```jl
julia> A = [ 2 5 4 2 3
            4 1 0.4 6 5
            0.9 4 8 7 10
            5 7 0.5 2 7
            2 6 8 9 4];

julia> B = [26; 24.4; 32.8; 30.5; 53];

julia> Xn = elimGaussNonPivoting(A,B)
5×1 Array{Float64,2}:
  2.0000000000000018
  2.999999999999999
  1.0000000000000013
  2.999999999999999
 -1.0000000000000009
```
return solusi SPL `X`
"""
function elimGaussNonPivoting(A,b)
    # Hitung ukuran matriks dan inisiasi solusi X
    n = length(b);
    X = zeros(n,1);
    # Buatlah matriks Gandeng [A|b]
    Aug = [A b];
    # Lakukan operasi baris dasar terhadap matriks gandeng
    for p = 1:n-1 # p menunjukkan kolom
        # Jika pivot bernilai nol, maka gagal.
        if Aug[p,p]==0
            error("Harus Pakai Pivoting")
        end
        # Lakukan eliminasi menggunakan operasi Eij(k)
        for i = p+1:n # i menunjukkan baris
            k = Aug[i,p]/Aug[p,p];
            Aug[i,:] = Aug[i,:] - k*Aug[p,:];
        end
    end
    # Pisahkan matriks gandeng, kemudian substitusi mundur.
    A = Aug[:,1:n];
    b = Aug[:,1+n];
    X = backsub(A,b);
end

"""
    elimGaussWithPivoting(A,B)

adalah fungsi untuk menyelesaikan SPL `Ax=B` dengan `A` adalah matriks sembarang
dan `B` adalah vektor ruas kanan menggunakan metode eliminasi Gauss dengan pivoting.

# Examples
```jl
julia> A = [ 2 5 4 2 3
            4 1 0.4 6 5
            0.9 4 8 7 10
            5 7 0.5 2 7
            2 6 8 9 4];

julia> B = [26; 24.4; 32.8; 30.5; 53];

julia> Xn = elimGaussWithPivoting(A,B)
5×1 Array{Float64,2}:
  1.9999999999999993
  3.000000000000001
  1.0
  3.0
 -1.0000000000000007
```
return solusi SPL `X`

"""
function elimGaussWithPivoting(A,b)
    # Hitung ukuran matriks dan inisiasi solusi X
    n = length(b);
    X = zeros(n,1);
    # Buatlah matriks Gandeng [A|b]
    Aug = [A b];
    # Lakukan operasi baris dasar terhadap matriks gandeng
    for p = 1:n-1
        # Lakukan operasi pindah baris untuk menentukan nilai pivot
        val,j = findmax(abs.(Aug[p:n,p]));
        C = Aug[p,:];
        Aug[p,:] = Aug[j+p-1,:];
        Aug[j+p-1,:] = C;
        # Jika pivot bernilai nol, maka gagal.
        if Aug[p,p] == 0;
            error("Matriks tidak punya solusi tunggal")
        end
        # Lakukan eliminasi menggunakan operasi Eij(k)
        for i = p+1:n
            k = Aug[i,p]/Aug[p,p];
            Aug[i,:] = Aug[i,:] - k*Aug[p,:];
        end
    end
    # Pisahkan matriks gandeng, kemudian substitusi mundur.
    A = Aug[:,1:n];
    b = Aug[:,1+n];
    X = backsub(A,b);
end

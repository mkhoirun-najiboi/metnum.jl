"""
    backsub(A,B)

adalah fungsi untuk menyelesaikan SPL `Ax=B` dengan `A` adalah matriks segitiga
atas dan `B` adalah vektor ruas kanan.

# Examples
```jldoctest
julia> A = [ 4 -1 2 3
             0 -2 7 -4
             0 0 6 5
             0 0 0 3];

julia> B = [20;-7; 4; 6];

julia> X = backsub(A,B)
4×1 Array{Float64,2}:
  3.0
 -4.0
 -1.0
  2.0
```
return solusi SPL `X`
"""
function backsub(A,B)
    # Hitung ukuran matriks dan inisiasi solusi X
    n = length(B);
    X = zeros(n,1);
    # Hitung nilai solusi X ke-n
    X[n] = B[n]/A[n,n];
    # Hitung nilai solusi X ke- n-1 sampai X ke-1
    for i = n-1:-1:1
        X[i] = (B[i] - A[i,i+1:n]'X[i+1:n])/A[i,i];
    end
    return X
end

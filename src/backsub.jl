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
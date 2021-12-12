# Metode Faktorisasi LU Dengan Pivoting 
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
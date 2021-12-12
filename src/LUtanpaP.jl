# Metode Faktorisasi LU Tanpa Pivoting
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
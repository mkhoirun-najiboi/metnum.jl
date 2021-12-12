# Metode Eliminasi Gauss Tanpa Pivoting
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
# Metode Eliminasi Gauss Dengan Pivoting
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
"""
    secant(f,p0,p1)

adalah fungsi untuk menyelesaikan `f(x)=0` menggunakan metode secant. Metode ini
memerlukan dua nilai tebakan awal `p0` dan `p1`.

# Examples
```jl
julia> f(x) = x^3-3*x+2;

julia> pk, flag, M = secant(f,-2.6,-2.4);

julia> pk
-2.0

julia> flag
0

julia> M
9×3 Array{Float64,2}:
 0.0  -2.6      NaN
 1.0  -2.4      NaN
 2.0  -2.1066     0.293401
 3.0  -2.02264    0.0839576
 4.0  -2.00151    0.0211303
 5.0  -2.00002    0.00148856
 6.0  -2.0        2.25138e-5
 7.0  -2.0        2.26855e-8
 8.0  -2.0        3.40616e-13
```
return solusi `p0` dengan `flag` bernilai 0 jika metode bisection berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, pk, error]`
"""
function secant(f,p0,p1)
    # Definisikan nilai toleransi, maksimum iterasi dan tebakan awal.
    delta = 10^-12;
    maxi = 100;
    M = [0 p0 NaN
         1 p1 NaN];
    p = [p0 p1];
    pk = p1;
    flag = 1;
    # Mulai langkah iterasi
    for k = 2:maxi
        # rumus metode secant
        pk=p[k]-f(p[k])*(p[k]-p[k-1])/(f(p[k])-f(p[k-1]));
        p = [p pk]
        # Hitung nilai galat mutlak dan relatif
        err=abs(p[k+1]-p[k]);
        rel=2*err/(abs(p[k]) + eps());
        M = [M; [k pk err]]
        # Kriteria penghentian iterasi jika galat memenuhi toleransi.
        if err<delta || rel<delta
            flag = 0;
            break
        end
    end
    return pk, flag, M
end

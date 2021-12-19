"""
    newtonRaphson(f,df,p0)

adalah fungsi untuk menyelesaikan `f(x)=0` menggunakan metode newton. Metode ini
memerlukan turunan dari fungsi `f(x)` yaitu `df` dan nilai tebakan awal `p0`.

# Examples
```jl
julia> f(x) = x^3-3*x+2;

julia> df(x) = 3*x^2-3;

julia> pk,flag,M = newtonRaphson(f,df,-2.4);

julia> pk
-2.0

julia> flag
0

julia> M
6×3 Array{Float64,2}:
 0.0  -2.4      NaN
 1.0  -2.07619    0.32381
 2.0  -2.0036     0.0725945
 3.0  -2.00001    0.00358742
 4.0  -2.0        8.58992e-6
 5.0  -2.0        4.91913e-11
```
return solusi `p0` dengan `flag` bernilai 0 jika metode bisection berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, pk, error]`
"""
function newtonRaphson(f,df,p0)
    # Definisikan nilai toleransi, maksimum iterasi dan tebakan awal
    delta = 10^-7;
    maxi = 100;
    pk = p0
    M = [0 pk NaN];
    # Mulai langkah iterasi
    flag = 1;
    for k = 2:maxi
        # Rumus metode Newton-Raphson
        pk1 = pk
        pk = pk1 - f(pk1)/df(pk1);
        # Hitung nilai galat mutlak dan relatif
        err = abs(pk-pk1);
        rel = 2*err/(abs(pk)+eps());
        M = [M; [k-1 pk err]]
        # Kriteria penghentian iterasi jika galat memenuhi toleransi.
        if err<delta || rel<delta
            flag = 0;
            break
        end
    end
    return pk, flag, M
end

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
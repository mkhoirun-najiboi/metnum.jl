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
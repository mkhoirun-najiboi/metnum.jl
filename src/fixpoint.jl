function fixpoint(g,p0)
    # Definisikan nilai toleransi, maksimum iterasi dan tebakan awal
    delta = 10^-7;
    maxi = 100;
    flag = 1;
    pn = p0
    M = [0 pn NaN];
    for n = 2:maxi
        pn1 = pn;
        pn = g(pn1);
        err = abs(pn-pn1);
        relerr = err/(abs(pn)+eps()); 
        M = [M; [n-1 pn err]]
        if (err<delta) || (relerr<delta)
            flag = 0; break
        end
    end
    return pn, flag, M
end
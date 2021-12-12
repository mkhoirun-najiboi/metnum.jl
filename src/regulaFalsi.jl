function regulaFalsi(f,a,b)
    delta = 10^-7; maxi = 100; flag = 1;
    M = Array{Float64}(undef, 0, 5);
    fa = f(a); fb = f(b);
    if fa*fb > 0
        c = "error : fa fb harus beda tanda";
        flag = 2; 
        return;
    end 
    for k = 1:maxi
        c = b-fb*(b-a)/(fb-fa);
        fc = f(c); 
        dx = min(c-a,b-c);
        M = [M ; [k-1 a c b fc] ]; 
        if fc == 0
            a = c;
            b = c;
        elseif fa*fc>0
            a = c;
            fa= fc;
        else
            b = c;
            fb= fc;
        end
        if abs(fc) < delta || abs(dx)< delta
            flag = 0; break;
        end
    end
    return c,flag,M
end
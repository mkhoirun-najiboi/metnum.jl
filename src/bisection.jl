function bisection(f,a,b)
    delta = 10^-7; maxi = 100; flag = 1;
    M = Array{Float64}(undef, 0, 5);
    fa = f(a); fb = f(b);
    if fa*fb>0
        c = "error: f(a) dan f(b) harus berbeda tanda";
        flag = 2; 
        return;
    end
    k = 1
    while k<=maxi
        c  = (b+a)/2; 
        fc = f(c); 
        M = [M; [k-1 a c b fc] ];
        if fc == 0
            a = c;
            b = c;
        elseif fa*fc>0
            a = c;
            fa = fc;
        else
            b = c;
            fb = fc;
        end
        if b-a < delta || abs(fc) < delta
            flag = 0;
            break;
        end
        k+=1
    end
    return c,flag, M
end
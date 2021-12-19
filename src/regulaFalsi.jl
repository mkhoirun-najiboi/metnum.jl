"""
    regulaFalsi(f,a,b)

adalah fungsi yang mencari nilai akar persamaan dari fungsi `f(x)` pada interval
di antara `a` dan `b`.

# Examples
```jldoctest
julia> f(x) = x*sin(x)-1;

julia> c, flag, M = regulaFalsi(f,0,2);

julia> c
1.1141571430336825

julia> flag
0

julia> M
4×5 Array{Float64,2}:
 0.0  0.0      1.09975  2.0      -0.0200192
 1.0  1.09975  1.12124  2.0       0.00983461
 2.0  1.09975  1.11416  1.12124   5.63036e-6
 3.0  1.09975  1.11416  1.11416   3.00226e-9
```
return solusi `c` dengan `flag` bernilai 0 jika metode regula falsi berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, a, c, b, f(c)]`

"""
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

"""
    bisection(f,a,b)

adalah fungsi yang mencari nilai akar persamaan dari fungsi `f(x)` pada interval
di antara `a` dan `b`.

# Examples
```jldoctest
julia> f(x) = x*sin(x)-1;

julia> c, flag, M = bisection(f,0,2);

julia> c
1.1141571998596191

julia> flag
0

julia> M
22×5 Array{Float64,2}:
  0.0  0.0      1.0      2.0      -0.158529
  1.0  1.0      1.5      2.0       0.496242
  2.0  1.0      1.25     1.5       0.186231
  3.0  1.0      1.125    1.25      0.015051
  ⋮
 18.0  1.11415  1.11415  1.11416  -3.22926e-6
 19.0  1.11415  1.11416  1.11416  -5.80313e-7
 20.0  1.11416  1.11416  1.11416   7.44159e-7
 21.0  1.11416  1.11416  1.11416   8.19227e-8
```
return solusi `c` dengan `flag` bernilai 0 jika metode bisection berhasil menemukan
solusi dan gagal jika tidak 0. Serta, matriks `M` yang berisi catatan proses tiap
iterasi `[k, a, c, b, f(c)]`

"""
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

"""
    rungekutta(f,a,b,y0,M)

berisi program untuk mencari solusi persamaan differensial `y'=f(t,y)`
dengan masalah nilai awal `y(a) = y0` pada interval `[a, b]`. Program ini secara default
berisi 5 masukan, yaitu fungsi `f(t,y)`, titik ujung interval penyelesaian `[a,b]`, nilai awal
`y0`, dan jumlah sub-interval `M`.

# Examples
```jldoctest
julia> rungekutta((t,y)->(t-y)/2,0,3,1,6)
7×2 Array{Float64,2}:
 0.0  1.0
 0.5  0.836426
 1.0  0.819628
 1.5  0.917142
 2.0  1.10368
 2.5  1.35956
 3.0  1.66943
```
return solusi masalah nilai awal `sol`.
"""
function rungekutta(f,a,b,y0,M)
    M = Int(M)
    h = (b-a)/M;
    T = a:h:b;
    Y = Array{Float64}(undef,length(T),1)
    Y[1] = y0;
    for k = 1:M
        f1 = f(T[k]     ,Y[k]        );
        f2 = f(T[k]+h/2 ,Y[k]+f1*h/2 );
        f3 = f(T[k]+h/2 ,Y[k]+f2*h/2 );
        f4 = f(T[k]+h   ,Y[k]+f3*h   );
        Y[k+1] = Y[k] + h/6*(f1+2*f2+2*f3+f4);
    end
    sol = [T Y];
end

"""
    euler(f,a,b,y0,M)
berisi program untuk mencari solusi persamaan differensial `y'=f(t,y)`
dengan masalah nilai awal `y(a) = y0` pada interval `[a, b]`. Program ini secara default
berisi 5 masukan, yaitu fungsi `f(t,y)`, titik ujung interval penyelesaian `[a,b]`, nilai awal
`y0`, dan jumlah sub-interval `M`.
# Example
```jl
julia> f(t,y) = (t-y)/2;

julia> sol = euler(f,0,3,1,6)
7×2 Array{Float64,2}:
 0.0  1.0
 0.5  0.75
 1.0  0.6875
 1.5  0.765625
 2.0  0.949219
 2.5  1.21191
 3.0  1.53394
```
return solusi masalah nilai awal `sol`.
"""
function euler(f,a,b,y0,M)
    M = Int(M)
    h = (b-a)/M;
    T = a:h:b;
    Y = Array{Float64}(undef,length(T),1)
    Y[1] = y0;
    # Mulai langkah iterasi euler
    for k = 1:M
        # Rumus iterasi euler
        Y[k+1] = Y[k]+h*f(T[k],Y[k])
    end
    sol = [T Y]
end

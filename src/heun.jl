"""
    heun(f,a,b,y0,M)

berisi program untuk mencari solusi persamaan differensial `y'=f(t,y)`
dengan masalah nilai awal `y(a) = y0` pada interval `[a, b]`. Program ini secara default
berisi 5 masukan, yaitu fungsi `f(t,y)`, titik ujung interval penyelesaian `[a,b]`, nilai awal
`y0`, dan jumlah sub-interval `M`.

# Examples
```jldoctest
julia> heun((t,y)->(t-y)/2,0,3,1,6)
7×2 Array{Float64,2}:
 0.0  1.0
 0.5  0.84375
 1.0  0.831055
 1.5  0.930511
 2.0  1.11759
 2.5  1.37311
 3.0  1.68212
```
return solusi masalah nilai awal `sol`.
"""
function heun(f,a,b,y0,M)
    M = Int(M)
    h = (b-a)/M;
    T = a:h:b;
    Y = Array{Float64}(undef,length(T),1)
    P = Array{Float64}(undef,length(T),1)
    Y[1] = y0;
    # Mulai langkah iterasi Heun
    for k = 1:M
        #% Rumus iterasi Heun
        P[k+1]=Y[k]+h*f(T[k],Y[k]);
        Y[k+1]=Y[k]+h/2*(f(T[k],Y[k])+f(T[k+1],P[k+1]));
    end
    sol = [T Y];
end

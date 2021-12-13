"""
    heun(f,a,b,y0,M)

heun adalah fungsi yang digunakan untuk menyelesaikan masalah nilai awal
dengan metode heun.
`f` adalah fungsi dari y'=f(t,y)
`a`, `b` adalah batas interval masalah
`y0` adalah nilai awal y(a)
`M` adalah banyaknya sub-interval

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

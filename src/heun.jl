"""
    heun(f,a,b,y0,M)

heun adalah fungsi yang digunakan untuk menyelesaikan masalah nilai awal
dengan metode heun.
`f` adalah fungsi dari y'=f(t,y)
`a`, `b` adalah batas interval masalah
`y0` adalah nilai awal y(a)
`M` adalah banyaknya sub-interval

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

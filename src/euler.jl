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
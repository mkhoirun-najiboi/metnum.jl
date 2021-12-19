"""
    rungekuttasistem(f,a,b,y0,M)
berisi program untuk mencari solusi sistem persamaan differensial `x'=f(t,z)`
dan `y'=g(t,z)` dengan masalah
nilai awal `x(a)=x0` dan `y(a)=y0` pada interval `[a,b]`. Program ini secara default
berisi 5 masukan, yaitu fungsi yang berisi `f(t,z)` dan `g(t,z)`, titik ujung interval `[a,b]`,
nilai awal `z0` yang berisi `[x0,y0]`, dan jumlah sub-interval `M`. Program dapat digunakan juga
untuk sistem persamaan diferensial lebih dari 2 persamaan.

# Example
```jl
julia> f(t,z) = [ z[1]+2*z[2] , 3*z[1]+2*z[2] ];

julia> sol = rungekuttasistem(f,0,0.2,[6,4],10)
11×3 Array{Float64,2}:
 0.0    6.0       4.0
 0.02   6.29355   4.53932
 0.04   6.61562   5.11949
 0.06   6.96853   5.74397
 ⋮
 0.14   8.74141   8.76532
 0.16   9.29021   9.6746
 0.18   9.88827  10.6561
 0.2   10.5396   11.7158
```
return solusi masalah nilai awal `sol`.
"""
function rungekuttasistem(f,a,b,y0,M)
  M = Int(M)
  h = (b-a)/M;
  T = a:h:b;
  Y = Array{Float64}(undef,M+1,length(y0))
  Y[1,:] = y0;
  for k = 1:M
    f1 = f(T[k]     ,Y[k,:]        );
    f2 = f(T[k]+h/2 ,Y[k,:]+f1*h/2 );
    f3 = f(T[k]+h/2 ,Y[k,:]+f2*h/2 );
    f4 = f(T[k]+h   ,Y[k,:]+f3*h   );
    Y[k+1,:] = Y[k,:] + h/6*(f1+2*f2+2*f3+f4);
  end
  sol = [T Y];
end

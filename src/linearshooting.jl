"""
    linearshooting(F1,F2,a,b,alpha,beta,M)
berisi program untuk mencari solusi persamaan differensial
`x''=f(t,x,x')` dengan masalah nilai batas `x(a)=alpha` dan `x(b)=beta` pada interval `[a,b]`.
Program ini secara default berisi 5 masukan, yaitu fungsi `F1` dan `F2` (sistem persamaan
differensial hasil transformasi dari `x''`), titik ujung interval `[a,b]`, nilai batas `[alpha,beta]`, dan
jumlah sub-interval `M`.
# Example
```jl
julia> F1(t,z) = [ z[2], 2*t/(1+t^2)*z[2]-2/(1+t^2)*z[1]+1 ];

julia> F2(t,z) = [ z[2], 2*t/(1+t^2)*z[2]-2/(1+t^2)*z[1] ];

julia> solusi = linearshooting(F1,F2,0,4,1.25,-0.95,40)
41×2 Array{Float64,2}:
 0.0   1.25
 0.1   1.29112
 0.2   1.31735
 0.3   1.32899
 ⋮
 3.7  -1.03334
 3.8  -1.01809
 3.9  -0.990472
 4.0  -0.95
```
return solusi masalah nilai batas `sol`.
"""
function linearshooting(F1,F2,a,b,alpha,beta,M)
  M = Int(M)
  Za = [alpha,0];
  sol = rungekuttasistem(F1,a,b,Za,M);
  U = sol[:,2];

  Za = [0,1];
  sol = rungekuttasistem(F2,a,b,Za,M);
  V = sol[:,2];

  T = sol[:,1];
  X = U + (beta-U[M+1])*V/V[M+1];
  solusi = [T X];
end

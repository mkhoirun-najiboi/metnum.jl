"""
    gelombang(f,g,a,b,c,m,n)
berisi program untuk mencari solusi persamaan differensial parsial yaitu persamaan gelombang.
Program ini secara
default berisi 7 masukan, yaitu fungsi nilai batas `f(x)` dan `g(x)`, batas solusi `a` dan `b`,
koefisien persamaan gelombang `c`, dan jumlah grid solusi x dan t yaitu `m` dan `n`.
# Example
```jl
julia> f(x) = sin(pi*x)+sin(2*pi*x);

julia> g(x) = 0;

julia> U = gelombang(f,g,1,0.5,2,11,11)
11×11 Adjoint{Float64,Array{Float64,2}}:
 0.0   0.896802      1.53884     1.76007   …  -0.14204      -0.363271   -0.278768   0.0
 0.0   0.769421      1.32844     1.53884       1.11022e-16  -0.210404   -0.181636   0.0
 0.0   0.431636      0.769421    0.948401      0.360616      0.181636    0.0683644  0.0
 0.0  -2.22045e-16   0.0515989   0.181636      0.769421      0.639384    0.363271   0.0
 ⋮                                         ⋱                                        ⋮
 0.0  -0.363271     -0.639384   -0.769421     -0.181636     -0.0515989   0.0        0.0
 0.0  -0.0683644    -0.181636   -0.360616     -0.948401     -0.769421   -0.431636   0.0
 0.0   0.181636      0.210404    0.0          -1.53884      -1.32844    -0.769421   0.0
 0.0   0.278768      0.363271    0.14204   …  -1.76007      -1.53884    -0.896802   0.0
```
return solusi persamaan gelombang `U`.
"""
function gelombang(f,g,a,b,c,m,n)
  h = a/(m-1);
  x = 0:h:1;
  k = b/(n-1);
  r = c*k/h;
  U = zeros(m,n);
  for i = 2:m-1
    U[i,1] = f(x[i]);
    U[i,2] = (1-r^2)*f(x[i]) + k*g(x[i]) + r^2/2*(f(x[i+1])+f(x[i-1]));
  end
  for j = 3:n
    for i = 2:(m-1)
      U[i,j] = (2-2*r^2)*U[i,j-1]+r^2*(U[i-1,j-1]+U[i+1,j-1])-U[i,j-2];
    end
  end
  U = U';
end

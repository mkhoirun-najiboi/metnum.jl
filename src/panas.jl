"""
    panas(f,c1,c2,a,b,c,m,n)
berisi program untuk mencari solusi
persamaan differensial parsial yaitu persamaan panas. Program ini secara default berisi
8 masukan, yaitu fungsi nilai awal `f(x)`, koefisien nilai batas `c1` dan `c2`, nilai batas solusi
`a` dan `b`, koefisien persamaan gelombang `c`, dan jumlah grid solusi x dan t yaitu `m` dan
`n`.
# Example
```jl
julia> f(x) = 4*x - 4*x.^2;

julia> U,R = panas(f,0,0,1,0.2,1,6,11);

julia> R
0.4999999999999999

julia> U
11×6 Adjoint{Float64,Array{Float64,2}}:
 0.0  0.64       0.96      0.96      0.64       0.0
 0.0  0.48       0.8       0.8       0.48       0.0
 0.0  0.4        0.64      0.64      0.4        0.0
 0.0  0.32       0.52      0.52      0.32       0.0
 ⋮                                              ⋮
 0.0  0.1375     0.2225    0.2225    0.1375     0.0
 0.0  0.11125    0.18      0.18      0.11125    0.0
 0.0  0.09       0.145625  0.145625  0.09       0.0
 0.0  0.0728125  0.117813  0.117813  0.0728125  0.0
```
return solusi persamaan panas `U` dan rasio `R`.
"""
function panas(f,c1,c2,a,b,c,m,n)
  h = a/(m-1);
  k = b/(n-1);
  r = c^2*k/h^2;
  U = zeros(m,n);
  U[1,:] .= c1;
  U[m,:] .= c2;
  U[2:m-1,1] = f.(h:h:(m-2)*h)';
  for j = 2:n
    for i = 2:m-1
      U[i,j]=(1-2*r)*U[i,j-1]+ r*(U[i-1,j-1]+U[i+1,j-1]);
    end
  end
  U=U';
  return U,r
end

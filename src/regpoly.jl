"""
    regpoly(X,Y,m)

adalah fungsi untuk mencari persamaan polinomial `y=C[m]x^m + ... + C[2]x^2 + C[1]x`
dari variabel dependen `Y` dan independen `X` dengan pangkat tertinggi `m`.

# Example
```jldoctest
julia> x = [-3, 0, 2, 4];

julia> y = [ 3, 1, 1, 3];

julia> C = regpoly(x,y,2)
3-element Array{Float64,1}:
  0.17846247712019525
 -0.19249542403904824
  0.8505186089078707
```
return nilai koefisien `C = [C[m],...,C[2],C[1]]`.
"""
function regpoly(X,Y,m)
  F = zeros(length(X),m+1)
  for k = 1: m+1
    F[:,k]=X.^(k-1);
  end
  A = F'*F;
  B = F'*Y;
  C = A\B;
  C = reverse(C,dims=1)
end

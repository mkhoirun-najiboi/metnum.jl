"""
    regpower(X,Y,m)

adalah fungsi untuk mencari persamaan pangkat `y=Ax^m` dari variabel dependen `Y` dan
independen `X` dengan pangkat `m`.

# Example
```jldoctest
julia> tk = [0.2,0.4,0.6,0.8,1.0];

julia> dk = [0.1960,0.7850,1.7665,3.1405,4.9075];

julia> A = regpower(tk,dk,2)
4.907303370786516
```
return nilai koefisien `A`.
"""
function regpower(X,Y,m)
  sumxy = (X.^m)'*Y
  sumx2 = (X.^m)'*(X.^m)
  A = sumxy/sumx2
end

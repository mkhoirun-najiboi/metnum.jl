"""
    rekursif(f,a,b,n)

berisi program untuk mencari integral numerik dari fungsi
`f(x)` pada interval `[a,b]` menggunakan Aturan Rekursif Trapesium.
Program ini secara default berisi 4 masukan, yaitu fungsi `f(x)`,
titik ujung interval `[a,b]` dan banyaknya iterasi `n`

# Example
```jl
julia> f(x) = 1/x;

julia> T = rekursif(f,1,5,5)
6×1 Array{Float64,2}:
 2.4
 1.8666666666666667
 1.6833333333333333
 1.628968253968254
 1.6144063238103485
 1.6106858960792332
```
return solusi integral numerik pada setiap iterasi `T`
"""
function rekursif(f,a,b,n)
  M = 1;
  h = b-a;
  T = Array{Float64}(undef,n+1,1)
  T[1] = h/2*(f(a)+f(b))
  for j = 1:n
    h = h/2;
    s = 0;
    for k = 1:M
      x = a+h*(2*k-1);
      s = s+f(x);
    end
    M = 2*M;
    T[j+1] = T[j]/2+h*s;
  end
  return T
end

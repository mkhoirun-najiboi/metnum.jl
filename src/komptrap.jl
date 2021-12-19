"""
    komptrap(f,a,b,M)
adalah fungsi yang digunakan untuk menghitung nilai integral `f` secara numerik
pada interval `[a,b]` menggunakan komposit trapesium dengan `M` sub-interval.

# Example
```jl
julia> a = 1; b = 6;

julia> f(x) = 2+sin(2*sqrt(x));

julia> y = komptrap(f,a,b,10)
8.193854565172531
```
return solusi `y`
"""
function komptrap(f,a,b,M)
  h = (b-a)/M;
  s = 0;
  for k = 1:M-1
    x = a+k*h;
    s = s+f(x);
  end
  y = h/2*(f(a)+f(b)+2*s)
end

"""
    kompsimp(f,a,b,M)
adalah fungsi yang digunakan untuk menghitung nilai integral `f` secara numerik
pada interval `[a,b]` menggunakan komposit Simpson dengan `M` sub-interval.

# Example
```jl
julia> a = 1; b = 6;

julia> f(x) = 2+sin(2*sqrt(x));

julia> y = kompsimp(f,a,b,10)
8.183447496636239
```
return solusi `y`
"""
function kompsimp(f,a,b,M)
  h = (b-a)/(2*M);
  s1 = 0;   s2 = 0;
  for k = 1:M
    x = a + (2*k-1)*h;
    s1 = s1 + f(x);
  end
  for k = 1:M-1
    x = a + (2*k)*h;
    s2 = s2 + f(x);
  end
  y = h/3*(f(a)+f(b)+4*s1+2*s2)
end

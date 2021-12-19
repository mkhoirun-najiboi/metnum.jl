"""
    romberg(f, a, b; maxi = 10, delta = 10^-9)
berisi program untuk mencari integral numerik dari fungsi `f(x)` pada
interval `[a,b]`. Program ini secara default berisi 3 masukan, yaitu fungsi `f(x)`, dan titik
ujung interval `[a,b]`. Secara default memiliki nilai maksimum iterasi `maxi=10` dan
toleransi error sebesar `delta = 1e-9`
# Examples
```jl
julia> f(x) = (x^2+x+1)*cos(x)
f (generic function with 1 method)

julia> sol, flag, err, R = romberg(f, 0, pi/2);

julia> sol
2.0381974270672245

julia> flag
0

julia> err
1.213065203842234e-10

julia> R
6×6 Array{Float64,2}:
 0.785398  0.0      0.0     0.0     0.0     0.0
 1.72681   2.04062  0.0     0.0     0.0     0.0
 1.96053   2.03844  2.0383  0.0     0.0     0.0
 2.01879   2.03821  2.0382  2.0382  0.0     0.0
 2.03335   2.0382   2.0382  2.0382  2.0382  0.0
 2.03698   2.0382   2.0382  2.0382  2.0382  2.0382

 julia> sol, flag, err, R = romberg(f, 0, pi/2, maxi=20, delta=1e-14);

 julia> sol
 2.038197427067236

 julia> flag
 0

 julia> err
 0.0

 julia> R
 8×8 Array{Float64,2}:
  0.785398  0.0      0.0     0.0     0.0     0.0     0.0     0.0
  1.72681   2.04062  0.0     0.0     0.0     0.0     0.0     0.0
  1.96053   2.03844  2.0383  0.0     0.0     0.0     0.0     0.0
  2.01879   2.03821  2.0382  2.0382  0.0     0.0     0.0     0.0
  2.03335   2.0382   2.0382  2.0382  2.0382  0.0     0.0     0.0
  2.03698   2.0382   2.0382  2.0382  2.0382  2.0382  0.0     0.0
  2.03789   2.0382   2.0382  2.0382  2.0382  2.0382  2.0382  0.0
  2.03812   2.0382   2.0382  2.0382  2.0382  2.0382  2.0382  2.0382
```
return nilai solusi integral numerik `sol`, status solusi
`flag`, hampiran galat `err`, dan matriks `R` yaitu matriks iterasi Romberg.
"""
function romberg(f, a, b; maxi = 10, delta = 10^-9)
  flag = 1;
  M = 1;
  h = b-a;
  err = 1;
  R = h/2*(f(a)+f(b));
  for J = 1:maxi
    # Rekursif Trapesium
    h = h/2;
    M = 2*M;
    s = 0;
    for p = 1:M/2
      x = a+h*(2*p-1);
      s = s+f(x);
    end
    R = [R zeros(size(R,1),1) ; R[J,1]/2 + h*s zeros(1,size(R,1))]
    # Aturan Romberg
    for K=1:J
      R[J+1,K+1]=R[J+1,K]+(R[J+1,K]-R[J,K])/(4^K-1);
    end
    err = abs.(R[J,J]-R[J+1,J+1]);
    if err<delta; flag=0; break; end
  end
  sol = R[end,end];
  return sol, flag, err, R
end

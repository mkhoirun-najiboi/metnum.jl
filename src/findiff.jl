"""
    findiff(p,q,r,a,b,alpha,beta,M)
berisi program untuk mencari solusi persamaan differensial
`x''=f(t,x,x')` dengan masalah nilai batas `x(a)=alpha` dan `x(b)=beta` pada interval `[a,b]`.
Program ini secara default berisi 7 masukan, yaitu fungsi `p(t)`, `q(t)` dan `r(t)` dari persamaan
`x''=p(t)x'+q(t)x+r(t)`, titik ujung interval `[a,b]`, nilai batas `[alpha,beta]`, dan
jumlah sub-interval `M`.
# Example
```jl
julia> p(t)= 2*t/(1+t^2);
       q(t)= -2/(1+t^2);
       r(t)= 1 + 0*t;

julia> solusi = findiff(p,q,r,0,4,1.25,-0.95,40)
41×2 Array{Float64,2}:
 0.0   1.25
 0.1   1.29077
 0.2   1.31665
 0.3   1.32791
 ⋮
 3.7  -1.03501
 3.8  -1.01924
 3.9  -0.991063
 4.0  -0.95
```
return solusi masalah nilai batas `sol`.
"""
function findiff(p,q,r,a,b,alpha,beta,M)
  h = (b-a)/M;
  T = a:h:b;
  T = T[2:end-1];
  #% Bangun Matriks B
  B = -h^2*r.(T);
  B[1] = B[1] + (1+h/2*p(T[1]))*alpha;
  B[end] = B[end] + (1-h/2*p(T[end]))*beta;
  # Bangun Matriks A - Bagian Diagonal
  Ad = 2 .+h^2*q.(T);
  #% Bangun Matriks A - Bagian Bawah Diagonal
  Tbawah = T[2:end];
  Abawah = -1 .-h/2*p.(Tbawah);
  #% Bangun Matriks A - Bagian Atas Diagonal
  Tatas = T[1:end-1];
  Aatas = -1 .+h/2*p.(Tatas);
  A = Tridiagonal(Abawah,Ad,Aatas)
  #% Selesaikan AX=B
  X = A\B;
  T = [a; T; b];
  X = [alpha; X ;beta];
  solusi = [T X];
end

"""
    lagrange(x, xd, yd)
adalah fungsi yang digunakan untuk mencari nilai interpolasi pada titik/vektor `x`, jika
diketahui suatu himpunan pasangan terurut `(xd,yd)`.

# Example
```jl
julia> xd = [1,2,3,5];

julia> yd = [1.06 1.12 1.34 1.78];

julia> lagrange(4,xd,yd)
1.6000000000000003

julia> lagrange([2.5,4,5.5],xd,yd)
3-element Array{Float64,1}:
 1.2175
 1.6000000000000003
 1.8025000000000002
```
return solusi hampiran interpolasi `y`.
"""
function lagrange(x, xd, yd)
  #% sedikit trik agar bisa menghitung nilai interpolasi pada banyak
  #% titik sekaligus. Intinya, tetap diproses satu titik demi titik.
  m=length(x);
  y=zeros(m)
  if m > 1
    for i=1:m
      # proses satu titik demi satu titik
      y[i]=lagrange(x[i], xd, yd);
    end
    return y
  end
  #% periksa jumlah titik dan tentukan derajat polinom
  ntitik = length(xd);
  n = ntitik-1;  #% derajat maksimum, sesuai jumlah titik yang ada
  #% hitung L_{n,k}(x) untuk k=1:(n+1)
  Ln=ones(1,ntitik);
  for i=1:ntitik
    for j=1:ntitik
      if i!=j
        Ln[i] = Ln[i] * (x-xd[j])/(xd[i]-xd[j]);;
      end
    end
  end
  #% hitung P_{N}(x)
  y = 0;
  for i=1:ntitik
    y = y + yd[i]*Ln[i];
  end
  return y
end

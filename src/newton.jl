"""
    newton(x, xd, yd)
adalah fungsi yang digunakan untuk mencari nilai interpolasi pada titik/vektor `x`, jika
diketahui suatu himpunan pasangan terurut `(xd,yd)`.

# Example
```jl
julia> xd = [1,2,3,5];

julia> yd = [1.06 1.12 1.34 1.78];

julia> y,D = newton(4,xd,yd);

julia> y
1.6

julia> D
4×4 Array{Float64,2}:
 1.06  0.0   0.0    0.0
 1.12  0.06  0.0    0.0
 1.34  0.22  0.08   0.0
 1.78  0.22  0.0   -0.02
```
return solusi hampiran interpolasi `y` dan matriks beda-terbagi `D`.
"""
function newton(x, xd, yd)
  m=length(x);
  y=zeros(m)
  if m > 1
    for i=1:m
      # proses satu titik demi satu titik
      y[i],D =newton(x[i], xd, yd);
    end
    return y,D
  end
  #% periksa jumlah titik dan tentukan derajat polinom
  ntitik = length(xd);
  #% hitung tabel beda-terbagi (divided-difference)
  D = zeros(ntitik,ntitik)
  D[:,1]=yd;          #% kolom pertama
  for j=2:ntitik      #% kolom ke-2 dan seterusnya
    for k=j:ntitik
        D[k,j] = (D[k,j-1]-D[k-1,j-1])/(xd[k]-xd[k-j+1]);
    end
  end
  #% hitung interpolasi Newton
  y = D[1,1]; s = 1;
  for i=2:ntitik
    s = s * (x-xd[i-1]);
    y = y + D[i,i]*s;
  end
  return y, D
end

"""
    bedaPusat(f, a; delta=10^-9)
adalah fungsi yang digunakan untuk mencari nilai turunan `f` pada titik `a`
menggunakan formula beda pusat. Secara default, digunakan nilai toleransi `delta=1e-9`

# Examples
```jl
julia> sol,flag,L = bedaPusat(sin,pi/3);

julia> sol
0.4999999999921733

julia> flag
0

julia> L
6×2 Array{Float64,2}:
 0.420735  NaN
 0.499167    0.0784316
 0.499992    0.000824583
 0.5         8.24996e-6
 0.5         8.24996e-8
 0.5         8.26006e-10

julia> sol,flag,L = bedaPusat(sin,pi/3,delta=1e-12);

julia> sol
0.4999999999588666

julia> flag
2

julia> L
8×2 Array{Float64,2}:
 0.420735  NaN
 0.499167    0.0784316
 0.499992    0.000824583
 0.5         8.24996e-6
 0.5         8.24996e-8
 0.5         8.26006e-10
 0.5         3.33067e-11
 0.5         3.33067e-10
```
return solusi `sol` dengan `flag` bernilai 0 jika toleransi terpenuhi, `flag`
bernilai 1 jika maksimum iterasi tercapai, `flag` bernilai 2 jika error minimum
telah tercapai namun tidak memenuhi toleransi. Serta, matriks `L` yang berisi
catatan tiap iterasi `[sol, error]`
"""
function bedaPusat(f, a; delta=10^-9)
  maxi = 15;
  flag = 1;
  h    = 1;
  D = (f(a+h)-f(a-h))/(2*h);
  E = NaN;
  sol = NaN
  for k = 1:maxi
    h = h/10;
    D = [D (f(a+h)-f(a-h))/(2*h)]
    E = [E abs(D[k+1]-D[k])]
    if E[k+1]<delta
      flag = 0;
      sol = D[end];
      break
    end
    if E[k+1]>E[k]
      sol  = D[k];
      flag = 2;
      break
    end
  end
  L = [D' E'];
  return sol, flag, L
end

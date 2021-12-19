"""
    richardson(f, a; delta = 1e-9)
adalah fungsi yang digunakan untuk mencari nilai turunan `f` pada titik `a`
menggunakan formula ekstrapolasi Richardson. Secara default, digunakan nilai
toleransi `delta=1e-9`

# Examples
```jl
julia> sol,flag,err,D = richardson(sin,pi/3);

julia> sol
0.49999999999998856

julia> flag
0

julia> err
3.323475383787411e-10

julia> D
5×5 Array{Float64,2}:
 0.420735  0.0       0.0       0.0  0.0
 0.479426  0.498989  0.0       0.0  0.0
 0.494808  0.499935  0.499998  0.0  0.0
 0.498699  0.499996  0.5       0.5  0.0
 0.499675  0.5       0.5       0.5  0.5

julia> sol,flag,err,D = richardson(sin,pi/3,delta=1e-12);

julia> sol
0.5000000000000013

julia> flag
0

julia> err
1.27675647831893e-14

julia> D
6×6 Array{Float64,2}:
 0.420735  0.0       0.0       0.0  0.0  0.0
 0.479426  0.498989  0.0       0.0  0.0  0.0
 0.494808  0.499935  0.499998  0.0  0.0  0.0
 0.498699  0.499996  0.5       0.5  0.0  0.0
 0.499675  0.5       0.5       0.5  0.5  0.0
 0.499919  0.5       0.5       0.5  0.5  0.5
```
return solusi `sol` dengan `flag` bernilai 0 jika toleransi terpenuhi, `flag`
bernilai 1 jika maksimum iterasi tercapai. Serta, error `err` dan matriks `D`
yang berisi tabel ekstrapolasi Richardson.
"""
function richardson(f, a; delta = 1e-9)
  maxi = 50;
  flag = 1;
  h = 1;
  D = (f(a+h)-f(a-h))/(2*h);
  err = NaN;
  for j=1:maxi
    h = h/2;
    D = [D zeros(size(D,1),1);
		(f(a.+h).-f(a.-h))./(2*h) zeros(1, size(D,1))];
    for k = 1:j
      D[j+1,k+1] = D[j+1,k] + (D[j+1,k]-D[j,k])/(4^k-1);
    end
    err = abs(D[j+1,j+1]-D[j,j]);
    if err<delta
      flag=0;
      break
    end
  end
  sol = D[end,end];
  return sol, flag, err, D
end

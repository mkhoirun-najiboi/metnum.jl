"""
    dirichlet(f1,f2,f3,f4,a,b,h)
berisi program untuk mencari solusi persamaan
differensial parsial yaitu persamaan Laplace dengan kondisi batas Dirichlet. Program
ini secara default berisi 8 masukan, yaitu fungsi nilai batas `f1`, `f2`, `f3`, dan
`f4`, nilai batas solusi `a` dan `b`, dan ukuran langkah `h`
# Example
```jl
julia> f1(x)= 20+0*x;

julia> f2(x)= 180+0*x;

julia> f3(y)= 80+0*y;

julia> f4(y)= 0*y;

julia> U = dirichlet(f1,f2,f3,f4,4,4,1)
5×5 Adjoint{Float64,Array{Float64,2}}:
  50.0   20.0      20.0      20.0     10.0
  80.0   55.7143   43.2143   27.1429   0.0
  80.0   79.6429   70.0      45.3571   0.0
  80.0  112.857   111.786    84.2857   0.0
 130.0  180.0     180.0     180.0     90.0
```
return solusi `U`.
"""
function dirichlet(f1,f2,f3,f4,a,b,h)
  maxi = 100;
  tol = 1e-7;
  n = Int(a/h)+1;
  m = Int(b/h)+1;
  U = (a.*(f1.(0).+f2.(0)).+b.*(f3.(0).+f4.(0)))./(2*a.+2*b).+ones(n,m);
  # Masalah nilai batas
  U[1:n,1]=f1.(0:h:(n-1)*h);
  U[1:n,m]=f2.(0:h:(n-1)*h);
  U[1,1:m]=f3.(0:h:(m-1)*h);
  U[n,1:m]=f4.(0:h:(m-1)*h);
  U[1,1]=(U[1,2]+U[2,1])/2;
  U[1,m]=(U[1,m-1]+U[2,m])/2;
  U[n,1]=(U[n-1,1]+U[n,2])/2;
  U[n,m]=(U[n-1,m]+U[n,m-1])/2;
  w = 4/(2+sqrt(4-(cos(pi/(n-1))+cos(pi/(m-1)))^2));
  err=1;   iter=0;
  while (err>tol)&&(iter<=maxi)
    err = 0;
    for j = 2:m-1
      for i = 2:n-1
        relx = w*(U[i,j+1]+U[i,j-1]+U[i+1,j]+U[i-1,j]-4*U[i,j])/4;
        U[i,j]=U[i,j]+relx;
        if err<=abs(relx)
          err=abs(relx);
        end
      end
    end
    iter=iter+1;
  end
  U = U';
end

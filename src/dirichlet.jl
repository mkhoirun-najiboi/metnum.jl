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
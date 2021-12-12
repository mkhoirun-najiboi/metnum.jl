function gelombang(f,g,a,b,c,m,n)
  h = a/(m-1);
  x = 0:h:1;
  k = b/(n-1);
  r = c*k/h;
  U = zeros(m,n);
  for i = 2:m-1
    U[i,1] = f(x[i]);
    U[i,2] = (1-r^2)*f(x[i]) + k*g(x[i]) + r^2/2*(f(x[i+1])+f(x[i-1]));
  end
  for j = 3:n
    for i = 2:(m-1)
      U[i,j] = (2-2*r^2)*U[i,j-1]+r^2*(U[i-1,j-1]+U[i+1,j-1])-U[i,j-2];
    end
  end
  U = U';
end
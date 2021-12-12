function panas(f,c1,c2,a,b,c,m,n)
  h = a/(m-1);
  k = b/(n-1);
  r = c^2*k/h^2;
  U = zeros(m,n);
  U[1,:] .= c1;
  U[m,:] .= c2;
  U[2:m-1,1] = f.(h:h:(m-2)*h)';
  for j = 2:n
    for i = 2:m-1
      U[i,j]=(1-2*r)*U[i,j-1]+ r*(U[i-1,j-1]+U[i+1,j-1]);
    end
  end
  U=U';
  return U,r
end
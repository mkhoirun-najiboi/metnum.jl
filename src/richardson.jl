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
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
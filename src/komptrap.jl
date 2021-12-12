function komptrap(f,a,b,M)
  h = (b-a)/M;
  s = 0; 
  for k = 1:M-1
    x = a+k*h;
    s = s+f(x);
  end 
  y = h/2*(f(a)+f(b)+2*s)
end
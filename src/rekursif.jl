function rekursif(f,a,b,n)
  M = 1;
  h = b-a;
  T = Array{Float64}(undef,n+1,1)
  T[1] = h/2*(f(a)+f(b))
  for j = 1:n
    h = h/2;
    s = 0;
    for k = 1:M
      x = a+h*(2*k-1);
      s = s+f(x);
    end
    M = 2*M;
    T[j+1] = T[j]/2+h*s;
  end
  return T
end
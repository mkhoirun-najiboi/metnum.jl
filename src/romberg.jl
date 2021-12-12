function romberg(f, a, b; maxi = 10, delta = 10^-9)
  flag = 1;
  M = 1;
  h = b-a;
  err = 1;
  R = h/2*(f(a)+f(b));
  for J = 1:maxi
    # Rekursif Trapesium
    h = h/2;
    M = 2*M;
    s = 0;
    for p = 1:M/2
      x = a+h*(2*p-1);
      s = s+f(x);
    end
    R = [R zeros(size(R,1),1) ; R[J,1]/2 + h*s zeros(1,size(R,1))]
    # Aturan Romberg
    for K=1:J
      R[J+1,K+1]=R[J+1,K]+(R[J+1,K]-R[J,K])/(4^K-1);
    end
    err = abs.(R[J,J]-R[J+1,J+1]);
    if err<delta; flag=0; break; end
  end
  sol = R[end,end];
  return sol, flag, err, R
end
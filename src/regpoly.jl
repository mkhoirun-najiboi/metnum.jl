function regpoly(X,Y,m)
  F = zeros(length(X),m+1)
  for k = 1: m+1
    F[:,k]=X.^(k-1);
  end
  A = F'*F;
  B = F'*Y;
  C = A\B;
  C = reverse(C,dims=1)
end

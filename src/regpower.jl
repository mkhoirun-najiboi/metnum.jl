function regpower(X,Y,m) 
  sumxy = (X.^m)'*Y
  sumx2 = (X.^m)'*(X.^m) 
  A = sumxy/sumx2
end
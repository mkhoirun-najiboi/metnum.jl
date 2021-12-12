function linearshooting(F1,F2,a,b,alpha,beta,M)
  M = Int(M)
  Za = [alpha,0];
  sol = rungekuttasistem(F1,a,b,Za,M);
  U = sol[:,2];
  
  Za = [0,1];
  sol = rungekuttasistem(F2,a,b,Za,M);
  V = sol[:,2];
  
  T = sol[:,1];
  X = U + (beta-U[M+1])*V/V[M+1];
  solusi = [T X];  
end
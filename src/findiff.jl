using LinearAlgebra
function findiff(p,q,r,a,b,alpha,beta,M)
  h = (b-a)/M;
  T = a:h:b;
  T = T[2:end-1];
  #% Bangun Matriks B
  B = -h^2*r.(T);
  B[1] = B[1] + (1+h/2*p(T[1]))*alpha;
  B[end] = B[end] + (1-h/2*p(T[end]))*beta;
  # Bangun Matriks A - Bagian Diagonal 
  Ad = 2 .+h^2*q.(T);
  #% Bangun Matriks A - Bagian Bawah Diagonal 
  Tbawah = T[2:end];
  Abawah = -1 .-h/2*p.(Tbawah);
  #% Bangun Matriks A - Bagian Atas Diagonal 
  Tatas = T[1:end-1];
  Aatas = -1 .+h/2*p.(Tatas);
  A = Tridiagonal(Abawah,Ad,Aatas)
  #% Selesaikan AX=B
  X = A\B;
  T = [a; T; b];
  X = [alpha; X ;beta];
  solusi = [T X];
end
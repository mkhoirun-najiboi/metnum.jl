function lagrange(x, xd, yd)
  #% sedikit trik agar bisa menghitung nilai interpolasi pada banyak 
  #% titik sekaligus. Intinya, tetap diproses satu titik demi titik.
  m=length(x);
  y=zeros(m)
  if m > 1
    for i=1:m
      # proses satu titik demi satu titik
      y[i]=lagrange(x[i], xd, yd);
    end
    return y
  end   
  #% periksa jumlah titik dan tentukan derajat polinom
  ntitik = length(xd);
  n = ntitik-1;  #% derajat maksimum, sesuai jumlah titik yang ada
  #% hitung L_{n,k}(x) untuk k=1:(n+1)
  Ln=ones(1,ntitik);
  for i=1:ntitik
    for j=1:ntitik
      if i!=j
        Ln[i] = Ln[i] * (x-xd[j])/(xd[i]-xd[j]);;
      end
    end
  end
  #% hitung P_{N}(x)
  y = 0;
  for i=1:ntitik
    y = y + yd[i]*Ln[i];
  end
  return y
end
function newton(x, xd, yd)
  m=length(x);
  y=zeros(m)
  if m > 1
    for i=1:m
      # proses satu titik demi satu titik
      y[i],D =newton(x[i], xd, yd);
    end
    return y,D
  end   
  #% periksa jumlah titik dan tentukan derajat polinom
  ntitik = length(xd);
  #% hitung tabel beda-terbagi (divided-difference)
  D = zeros(ntitik,ntitik)
  D[:,1]=yd;          #% kolom pertama
  for j=2:ntitik      #% kolom ke-2 dan seterusnya
    for k=j:ntitik
        D[k,j] = (D[k,j-1]-D[k-1,j-1])/(xd[k]-xd[k-j+1]);
    end
  end 
  #% hitung interpolasi Newton
  y = D[1,1]; s = 1;
  for i=2:ntitik
    s = s * (x-xd[i-1]);
    y = y + D[i,i]*s;
  end
  return y, D
end
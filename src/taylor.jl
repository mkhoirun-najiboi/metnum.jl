function taylor(df,a,b,y0,M)
  M = Int(M)
  h = (b-a)/M;
  T = a:h:b;
  Y = Array{Float64}(undef,length(T),1) 
  Y[1]=y0;
  for j = 1:M
    D = df(T[j],Y[j]);
    Y[j+1]=Y[j]+h*(D[1]+h*(D[2]/2+h*(D[3]/6+h*D[4]/24)));
  end
  sol = [T Y];
end
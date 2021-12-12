function rungekuttasistem(f,a,b,y0,M)
  M = Int(M)
  h = (b-a)/M;
  T = a:h:b;
  Y = Array{Float64}(undef,M+1,length(y0))
  Y[1,:] = y0;
  for k = 1:M
    f1 = f(T[k]     ,Y[k,:]        );
    f2 = f(T[k]+h/2 ,Y[k,:]+f1*h/2 );
    f3 = f(T[k]+h/2 ,Y[k,:]+f2*h/2 );
    f4 = f(T[k]+h   ,Y[k,:]+f3*h   );
    Y[k+1,:] = Y[k,:] + h/6*(f1+2*f2+2*f3+f4);
  end
  sol = [T Y];
end
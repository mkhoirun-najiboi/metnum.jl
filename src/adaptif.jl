function adaptif(f, a, b; delta=10^-9)
  iterating = 0;
  done = 1;
  SRvec = srule(f,a,b,delta);
  SRmat = SRvec;
  m = 1;
  state = iterating;
  while(state==iterating)
    n = m;
    for j = n:-1:1
      p = j;
      SR0vec = SRmat[p,:];
      err = SR0vec[5];
      delta = SR0vec[6];
      if (delta<=err)
        state=done;
        SR1vec=SR0vec;
        SR2vec=SR0vec;
        a = SR0vec[1];
        b = SR0vec[2];
        c =(a+b)/2;
        err=SR0vec[5];
        delta=SR0vec[6];
        delta2=delta/2;
        SR1vec=srule(f,a,c,delta2);
        SR2vec=srule(f,c,b,delta2);
        err = abs(SR0vec[3]-SR1vec[3]-SR2vec[3])/10;
        if err<delta
          SRmat[p,:]=SR0vec;
          SRmat[p,4]=SR1vec[3]+SR2vec[3];
          SRmat[p,5]=err;
        else
          SRmat = [SRmat; zeros(1,6)]
          SRmat[p+1:m+1,:] = SRmat[p:m,:];
          m=m+1;
          SRmat[p,:]=SR1vec;
          SRmat[p+1,:]=SR2vec;
          state=iterating;
        end
      end
    end
  end
  sol = sum(SRmat[:,4]);
  err = sum(abs.(SRmat[:,5]));
  SRmat = SRmat[1:m,:];
  return sol, err, SRmat
end
                    
function srule(f,a0,b0,delta0)
  h = (b0-a0)/2;
  C = zeros(1,3)
  C = f.([a0 (a0+b0)/2 b0]);
  S = h*(C[1]+4*C[2]+C[3])/3;
  S2= S;
  delta1=delta0;
  err=delta0;
  Z = [a0 b0 S S2 err delta1];
  return Z
end
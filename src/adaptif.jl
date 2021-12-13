"""
adaptif(f, a, b; delta=10^-9)

adalah fungsi yang digunakan untuk mencari nilai integral dari fungsi `f`
pada interval `[a,b]`. Secara default toleransi yang digunakan adalah `delta=1e-9`

# Examples
```jldoctest

julia> sol,err,SRmat = adaptif(f,0,4,delta=0.00001);

julia> sol
-1.5487882341253174

julia> err
2.9680861581356417e-6

julia> SRmat
20×6 Array{Float64,2}:
 0.0     0.0625   0.0228717  …  1.52153e-8  1.5625e-7
 0.0625  0.125    0.0594867     1.31649e-8  1.5625e-7
 0.125   0.1875   0.084342      1.13722e-8  1.5625e-7
 0.1875  0.25     0.0996986     9.80639e-9  1.5625e-7
 0.25    0.375    0.216719      2.50545e-7  3.125e-7
 0.375   0.5      0.206462   …  1.84024e-7  3.125e-7
 0.5     0.625    0.171505      1.33805e-7  3.125e-7
 0.625   0.75     0.124333      9.61076e-8  3.125e-7
 0.75    0.875    0.0732445     6.79932e-8  3.125e-7
 0.875   1.0      0.0235284     4.71839e-8  3.125e-7
 1.0     1.125   -0.0216607  …  3.19179e-8  3.125e-7
 1.125   1.25    -0.060651      2.0837e-8   3.125e-7
 1.25    1.5     -0.210811      3.17143e-7  6.25e-7
 1.5     2.0     -0.60551       3.19486e-8  1.25e-6
 2.0     2.25    -0.319856      8.10596e-8  6.25e-7
 2.25    2.5     -0.300617   …  8.30109e-8  6.25e-7
 2.5     2.75    -0.270099      7.07089e-8  6.25e-7
 2.75    3.0     -0.234747      5.4474e-8   6.25e-7
 3.0     3.5     -0.363888      1.03699e-6  1.25e-6
 3.5     4.0     -0.243134      4.10779e-7  1.25e-6

"""

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

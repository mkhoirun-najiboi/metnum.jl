"""
    adaptif(f, a, b; delta=10^-9)

adalah fungsi yang digunakan untuk mencari nilai integral dari fungsi `f`
pada interval `[a,b]`. Secara default toleransi yang digunakan adalah `delta=1e-9`

# Examples
```jldoctest
julia> f(x) = 13*(x-x.^2).*exp(-3*x/2);

julia> sol,err,SRmat = adaptif(f,0,4,delta=0.00001);

julia> sol
-1.5487882341253174

julia> err
2.9680861581356417e-6

julia> SRmat 
20×5 Array{Float64,2}:
 0.0     0.0625   0.0228718  1.52153e-8  1.5625e-7
 0.0625  0.125    0.0594869  1.31649e-8  1.5625e-7
 0.125   0.1875   0.0843421  1.13722e-8  1.5625e-7
 0.1875  0.25     0.0996987  9.80639e-9  1.5625e-7
 ⋮
 2.5     2.75    -0.2701     7.07089e-8  6.25e-7
 2.75    3.0     -0.234747   5.4474e-8   6.25e-7
 3.0     3.5     -0.363898   1.03699e-6  1.25e-6
 3.5     4.0     -0.243138   4.10779e-7  1.25e-6
```
`SRmat` adalah matriks yang berisi sub-interval (kolom 1 dan 2), nilai integral
pada sub-interval (kolom 3), galat integral numerik (kolom 4), dan toleransi pada
sub-interval (kolom 5).

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
  SRmat = SRmat[1:m,[1,2,4,5,6]];
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

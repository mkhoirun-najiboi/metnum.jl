"""
    rkf45(f,a,b,y0,M,delta)
berisi program untuk mencari solusi persamaan
differensial `y' = f(t,y)` dengan `y(a) = y0` pada interval `[a, b]`. Program ini secara
default berisi 6 masukan, yaitu fungsi `f(t,y)`,
titik ujung interval penyelesaian `[a,b]`,
nilai awal `y0`, jumlah sub-interval `M` dan nilai toleransi error `delta`

# Example
```jl
julia> f(t,y) = 1+y^2;

julia> sol = rkf45(f,0,1.4,0,7,2e-5)
14×2 Array{Float64,2}:
 0.0     0.0
 0.2     0.20271
 0.6     0.684165
 0.8     1.02968
 ⋮
 1.35    4.45586
 1.375   5.04273
 1.3875  5.39534
 1.4     5.79895
```
return solusi masalah nilai awal `sol`.
"""
function rkf45(f,a,b,y0,M,delta)
  M = round(M)
  a2=1/4;b2=1/4;a3=3/8;b3=3/32;c3=9/32;
  a4=12/13;b4=1932/2197;c4=-7200/2197;d4=7296/2197;
  a5=1;b5=439/216;c5=-8;d5=3680/513;e5=-845/4104;
  a6=1/2;b6=-8/27;c6=2;d6=-3544/2565;e6=1859/4104;f6=-11/40;
  r1=1/360;r3=-128/4275;r4=-2197/75240;r5=1/50;r6=2/55;
  n1=25/216;n3=1408/2565;n4=2197/4104;n5=-1/5;
  big=1e15;
  h=(b-a)/M;
  hmin=h/64;
  hmax=h*64;
  maxi=200;
  j=1;
  Y = y0;
  T = a;
  br= b-0.001*abs(b);
  err=NaN
  while T[j]<b
    if (T[j]+h)>br;h=b-T[j];end
    #% Hitung koefisien
    k1=h*f(T[j],Y[j]);
    y2=Y[j]+b2*k1;
    k2=h*f(T[j]+a2*h,y2);
    y3=Y[j]+b3*k1+c3*k2;
    k3=h*f(T[j]+a3*h,y3);
    y4=Y[j]+b4*k1+c4*k2+d4*k3;
    k4=h*f(T[j]+a4*h,y4);
    y5=Y[j]+b5*k1+c5*k2+d5*k3+e5*k4;
    k5=h*f(T[j]+a5*h,y5);
    y6=Y[j]+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5;
    k6=h*f(T[j]+a6*h,y6);
    err=abs(r1*k1+r3*k3+r4*k4+r5*k5+r6*k6);
    ynew=Y[j]+n1*k1+n3*k3+n4*k4+n5*k5;
    #% Perbarui ukuran langkah
    if (err<delta) || (h<2*hmin)
      Y = [Y; ynew];
      if (T[j]+h)>br
        T = [T; b];
      else
        T = [T; T[j]+h];
      end
      j=j+1;
    end
    if (err==0)
      s=0;
    else
      s=0.84*(delta*h/err)^(0.25);
    end
    if (s<0.75)&&(h>2*hmin);h = h/2;end
    if (s>1.50)&&(2*h<hmax);h = 2*h;end
    if abs(Y[j])>big || maxi==j;break;end
  end
  sol = [T Y];
  return sol
end

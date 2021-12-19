"""
    taylor(df,a,b,y0,M)

berisi program untuk mencari solusi persamaan differensial `y'=f(t,y)`
dengan masalah nilai awal `y(a) = y0` pada interval `[a, b]`. Program ini secara default
berisi 5 masukan, yaitu fungsi `df` yang berisi nilai dari fungsi y', y'', y''', y'''',
titik ujung interval penyelesaian `[a,b]`, nilai awal
`y0`, dan jumlah sub-interval `M`.
# Examples
```jldoctest
julia> function df(t,y)
         y1 = (t-y)/2;
         y2 = (2-t+y)/4;
         y3 = (-2+t-y)/8;
         y4 = (2-t+y)/16;
         z = [y1,y2,y3,y4];
       end

julia> sol1 = taylor(df,0,3,1,6)
7×2 Array{Float64,2}:
 0.0  1.0
 0.5  0.836426
 1.0  0.819628
 1.5  0.917142
 2.0  1.10368
 2.5  1.35956
 3.0  1.66943
```
return solusi masalah nilai awal `sol`.
"""
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

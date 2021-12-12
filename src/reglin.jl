"""
    reglin(X,Y)

regresi linear untuk variabel dependen `Y` dan independen `X` dengan persamaan
    y=Ax+B
return nilai `A` dan `B`

#Example
```jldoctest
julia> x = [-1, 0, 1, 2, 3, 4, 5, 6];
julia> y = [10, 9, 7, 5, 4, 3, 0, -1];
julia> A,B = reglin(x,y)
(-1.6071428571428572, 8.642857142857142)
"""
function reglin(X,Y)
    xmean = mean(X);
    ymean = mean(Y);
    # Hitung nilai jumlah dari xy dan x^2
    sumxy = (X.-xmean)'*(Y.-ymean)
    sumx2 = (X.-xmean)'*(X.-xmean)
    # Hitung nilau koefisien garis regresi linear Y=Ax+B
    A = sumxy/sumx2;
    B = ymean .- A*xmean;
    return A,B
end

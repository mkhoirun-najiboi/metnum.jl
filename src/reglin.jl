"""
    reglin(X,Y)

regresi linear untuk variabel dependen Y dan independen X

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

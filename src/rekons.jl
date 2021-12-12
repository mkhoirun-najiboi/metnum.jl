using LinearAlgebra
function rekons(A,B,Xawal::Array{Int64,1})
    delta = 10^-7;
    maxi = 100;
    flag = 1; 
    X = Xawal;
    Xlama = X;
    M = [0 X' NaN]; 
    for k = 1:maxi
        for i = 1:length(B)
            Xlama = X; 
            X = Xlama + A[i,:]*(B[i]-A[i,:]'Xlama)/(A[i,:]'A[i,:])
        end
        err = norm(X-Xlama) 
        M = [M; [k X' err] ];
        if err<delta
            flag = 0;
            break
        end
    end
    err = M[end,end]
    return X, flag, err, M
end
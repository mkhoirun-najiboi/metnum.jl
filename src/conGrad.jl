using LinearAlgebra
function conGrad(A,B,Xawal::Array{Int64,1})
    if ~(isposdef(A) && issymmetric(A))
        error("matriks A harus simetrik definit positif")
    end
    delta = 10^-7;
    maxi = 100;
    flag = 1; 
    X = Xawal;
    r = B-A*X
    d = r
    M = [0 X' NaN]; 
    for k = 1:maxi
        Xlama = X
        rlama = r
        dlama = d
        X = Xlama + ((rlama'rlama)/(dlama'A*dlama))*dlama
        r = rlama - ((rlama'rlama)/(dlama'A*dlama))*A*dlama
        d = r - ((r'r)/(rlama'rlama))*dlama
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
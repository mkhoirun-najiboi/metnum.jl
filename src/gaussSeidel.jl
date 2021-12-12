using LinearAlgebra
function gaussSeidel(A,B,Xawal::Array{Int64,1})  
    delta = 10^-7;
    maxi = 100;
    flag = 1; 
    D = Diagonal(diag(A))
    R = tril(A);
    U = triu(A,1); 
    X = Xawal; 
    M = [0 X' NaN]; 
    for k = 1:maxi
        Xlama = X;
        X = inv(R)*(B - U*Xlama); 
        err = norm(X-Xlama);
        M = [M; [k X' err] ];
        if err<delta
            flag = 0;
            break
        end
    end
    err = M[end,end]
    return X, flag, err, M
end
function [SigmaZ, svp] = ClosedNNFN(SigmaY, lambda, alpha, C)
norm_inf = SigmaY(1);
Lambda = lambda*2*sqrt(C);

if (norm_inf>Lambda)
    z = max(abs(SigmaY)-Lambda, 0);

    norm_2 = norm(z, 2);
    K = 1 + alpha*Lambda/norm_2;
    SigmaZ = K*z + sqrt(C);
   
    over = find(SigmaY>(K*Lambda-sqrt(C))/(K-1));
    SigmaZ(over) = SigmaY(over) - eps;
    ind = find(z > 0); 
    svp = length(ind);
    SigmaZ = SigmaZ(ind);
    
else 
    SigmaZ = [];    svp = 0;
end
return;

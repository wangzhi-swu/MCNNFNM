function [ Y_hat, W_hat] = NNFN_Estimation( NL_mat, Sigma_arr, CurPat, Par )

    Y_hat = zeros(size(CurPat));
    W_hat = zeros(size(CurPat));
    for i   =  1 : length(Par.SelfIndex) % For each keypatch group
        Y   =   CurPat(:, NL_mat(1:Par.nlsp,i)); % Non-local similar patches to the keypatch
        mY  =   repmat(mean( Y, 2 ),1,Par.nlsp);
        Y   =   Y-mY;
        
        X 	=   NNFN_ADMM( Y, Sigma_arr(:, Par.SelfIndex(i)), Par); 

        Y_hat(:,NL_mat(1:Par.nlsp,i))  = Y_hat(:,NL_mat(1:Par.nlsp,i)) + X+mY;
        W_hat(:,NL_mat(1:Par.nlsp,i))  = W_hat(:,NL_mat(1:Par.nlsp,i)) + ones(Par.ps2ch, Par.nlsp);
    end
end

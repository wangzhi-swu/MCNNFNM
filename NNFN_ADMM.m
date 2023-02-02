function  [Z] =  NNFN_ADMM( Y, NSig, Par )

if ~isfield(Par, 'maxIter')
    Par.maxIter = 10;
end
if ~isfield(Par, 'rho')
    Par.rho = 1;
end
if ~isfield(Par, 'mu')
    Par.mu = 1;
end
if ~isfield(Par, 'display')
    Par.display = true;
end
% Initializing optimization variables
% Intialize the weight matrix W
mNSig = min(NSig);

W = (mNSig+eps) ./ (NSig+eps);
% Initializing optimization variables
X = zeros(size(Y));
Z = zeros(size(Y));
A = zeros(size(Y));
%% Start main loop
iter = 0;
% T1 = zeros(Par.maxIter,1);T2 = zeros(Par.maxIter,1); T3 = zeros(Par.maxIter,1);
while (iter < Par.maxIter)
    iter = iter + 1;
    
    % update X, fix Z and A
    % min_{X} ||W * Y - W * X||_F^2 + 0.5 * rho * ||X - Z + 1/rho * A||_F^2
    X = diag(1 ./ (W.^2 + 0.5 * Par.rho)) * (diag(W.^2) * Y + 0.5 * Par.rho * Z - 0.5 * A);
    
    % update Z, fix X and A
    [U, SigmaY, V] = svd(full(X + A/Par.rho), 'econ');
    [SigmaZ , svp] = ClosedNNFN(diag(SigmaY), Par.lambda, Par.alpha, 2*Par.Constant*sqrt(Par.nlsp)*mNSig^2/Par.rho);
    Z = U(:, 1:svp) * diag(SigmaZ) * V(:, 1:svp)';
    
    % update the multiplier A, fix Z and X
    A = A + Par.rho*(X - Z);
    
    Par.rho = min(1e4, Par.mu * Par.rho);
end
return;

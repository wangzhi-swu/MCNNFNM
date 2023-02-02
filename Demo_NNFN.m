clear;
Original_image_dir  =   'kodak_color/'; 
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

nSig = [30 10 50];

for Search_Win = [20]
P_S = floor(Search_Win/3);
for LAMBDA= [0.86]
for RHO = [0.86]
for ALPHA=[1.9]
        
% parameters for denoising
Par.nSig      =   nSig;          % STD of the noise image
Par.win       =   Search_Win;    % Non-local patch searching window
Par.Constant  =   2 * sqrt(2);   % Constant num for the weight vector
Par.Innerloop =   2;             % InnerLoop Num of between re-blockmatching
Par.ps        =   P_S;           % Patch size, larger values would get better performance, but will be slower
Par.step      =   P_S-1;         % The step between neighbor patches, smaller values would get better performance, but will be slower
Par.Iter      =   10;            % total iter numbers
Par.display   =   true;

% parameters for ADMM
Par.maxIter   =   10;
Par.delta     =   0.1;             % iterative regularization parameter
Par.mu        =   1.002;
Par.rho       =   RHO;             % In final version, this parameter will be changed

    % this parameter is not finally determined yet
Par.lambda    =   LAMBDA;          % for different noise levels, this parameter should be tuned to achieve better performance


% record all the results in each iteration
Par.PSNR = zeros(Par.Iter, im_num, 'single');
for i = 1:im_num


    Par.alpha = ALPHA;
    Par.mu = Par.mu;
    Par.image   =   i;
    Par.nSig0   =   nSig;
    Par.nlsp    =   70;
    Par.I = double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
    S = regexp(im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);
   
    Par.nim = zeros([h, w, ch]); % add noise
    for c = 1:ch
        randn('seed',0);
        Par.nim(:, :, c) = Par.I(:, :, c) + Par.nSig0(c) * randn(size(Par.I(:, :, c)));
    end

    fprintf('%s  lambda=%.2f  rho=%.2f  mu=%.4f  alpha=%.2f  Win=%d:\n',im_dir(i).name, Par.lambda, Par.rho, Par.mu, Par.alpha, Par.win);
    PSNR  =   csnr( Par.nim, Par.I, 0, 0 );
    fprintf('The initial value of PSNR = %f\n', PSNR);
    
    time0 = clock;
    [im_out, Par] = NNFN_Denoising( Par.nim, Par.I, Par );
    im_out(im_out>255) = 255;
    im_out(im_out<0.0) = 0.0;
    fprintf('Total elapsed time = %f s\n',(etime(clock, time0)));
%     imwrite(im_out./255, fullfile([im_dir(i).name(6:7), '_lambda_', num2str(Par.lambda), '.png']));


end
end
end
end
end


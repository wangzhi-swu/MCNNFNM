function [rI, Par]   =  NNFN_Denoising( nI, I, Par )
rI  =  nI;   % Estimated Image
[h, w, ch]  = size(rI);
Par.h = h;
Par.w = w;
Par.ch = ch;

Par = SearchNeighborIndex( Par ); 
                    % 主要是 向Par增加了NumIndex、NeighborIndex和SelfIndex
% sig_Y = zeros(Par.Iter, 1);
% noisy image to patch
NoiPat   =	Image2Patch( nI, Par ); % NoiPat(:,i) := 第i块patch的108个位置的灰度值(108=3*6*6, row first)
Par.TolN =  size(NoiPat, 2); % patch 总数
Sigma_arrCh = zeros(Par.ch, Par.TolN); % 
for iter = 1 : Par.Iter
    Par.iter = iter;
    % iterative regularization 
    rI =	rI + Par.delta * (nI - rI); % 缝针
    % image to patch
    CurPat =	Image2Patch( rI, Par );
    % estimate local noise variance
    for c = 1:Par.ch
%         if(iter == 1)
%             TempSigma_arrCh = sqrt(max(0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
%             %             TempSigma_arrCh = sqrt(abs(repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
%         else
%             TempSigma_arrCh = Par.lambda*sqrt(max(0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
%             %             TempSigma_arrCh = Par.lambda*sqrt(abs(repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
%         end
        TempSigma_arrCh = sqrt(max(0, repmat(Par.nSig0(c)^2, 1, size(CurPat, 2)) - mean((NoiPat((c-1)*Par.ps2+1:c*Par.ps2, :) - CurPat((c-1)*Par.ps2+1:c*Par.ps2, :)).^2)));
        Sigma_arrCh((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat(TempSigma_arrCh, [Par.ps2, 1]);
    end
%     sig_Y(iter) = mean(Sigma_arrCh([1,Par.ps2+1,2*Par.ps2+1],1));
    
    if (mod(iter-1, Par.Innerloop) == 0) % iter = 1, 3, 5, ...
        Par.nlsp = Par.nlsp - 10;  % Lower Noise level, less NL patches
        NL_mat  =  Block_Matching(CurPat, Par); % Caculate Non-local similar patches for each
        if iter>2 % iter=3, 5, 7, ...
            Par.mu = max(1+eps, Par.mu-0.001);
        end
    end
   
    [Y_hat, W_hat]  =  NNFN_Estimation( NL_mat, Sigma_arrCh, CurPat, Par );   % Estimate all the patches
    
    rI = PGs2Image(Y_hat, W_hat, Par);
    rI(rI>255) = 255;
    rI(rI<0.0) = 0.0;
    PSNR   =  csnr( I, rI, 0, 0 );
%     SSIM   =  cal_ssim( I, rI, 0, 0 );
    fprintf( 'Iter=%2.0f, PSNR = %2f\n', iter, PSNR );
    Par.PSNR(iter, Par.image)  =  PSNR;
%     Par.SSIM(iter, Par.image)  =  SSIM;
    
%     Par.lambda = Par.lambda/log(Par.rho); % 减少lambda
    
    if (iter>1 && Par.PSNR(iter-1, Par.image)>PSNR)
%         plot(1:iter, sig_Y(1:iter), 'p');
        break;
    end
end
return;






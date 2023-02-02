function       X = Image2Patch( im_out, par )
% record the non-local patch set and the index of each patch in
% of seed patches in image
im_out     =  single(im_out); % double -> single ( double占8bits, single占4。所以用single会更快一些 )
X          =  zeros(par.ps2ch, par.maxrc, 'double');
k          =  0;
for channel = 1:par.ch
    for i = 1:par.ps
        for j = 1:par.ps
            k    =  k+1;
            blk  =  im_out(i:end-par.ps+i,j:end-par.ps+j, channel);
            X(k,:) = blk(:)';
        end
    end
end

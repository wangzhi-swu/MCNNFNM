function       blk_arr = Block_Matching( X, par)
blk_arr   =  zeros(par.nlsp, par.lenrc, 'single');

for  i  =  1 : par.lenrc
    seed = X(:, par.SelfIndex(i));
    neighbor = X(:, par.NeighborIndex(1:par.NumIndex(i), i));
    dis = sum(bsxfun(@minus, neighbor, seed).^2, 1);
    [~,ind]   =  sort(dis);
    indc = par.NeighborIndex( ind( 1:par.nlsp ), i );
    blk_arr(:, i) = indc;
end
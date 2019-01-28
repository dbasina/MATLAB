function pts_out=fracSingleM(pts_tar,mat_p1,iter,PERM)
nT      = size(pts_tar,1);
D       = max(size(mat_p1));
cellOut = cell(nT,1);
    parfor i=1:nT
        if(PERM)
            mat_p1T = reshape(mat_p1(randperm(numel(mat_p1))),size(mat_p1));
        else
            mat_p1T = mat_p1;    
        end
        IDm     = find(mat_p1T);
        tmpSize = numel(IDm);
        [x y z] = ind2sub(size(mat_p1T),IDm);
        pts_tmp = [(([x y z] -0.5 - D/2)/D)/(D^(iter-1)) mat_p1T(IDm) repmat(iter,tmpSize,1)];
        pts_tmp(:,1:3)  = pts_tmp(:,1:3) + repmat(pts_tar(i,1:3),tmpSize,1);
        pts_tmp(:,4)    = pts_tmp(:,4)  .* repmat(pts_tar(i,4),tmpSize,1);
        cellOut{i} = pts_tmp;
    end
pts_out=cell2mat(cellOut);
end
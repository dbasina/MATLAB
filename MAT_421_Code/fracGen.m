function pts_out=fracGen(varargin)
% Generates a multi/monofractal point distribution based on the user supplied 
% or build-in template matrix. 
%
% syntax 1: fracGen(fractalName)
%fractalName    = name of the build-in fractal. An invalid name will trigger 
%the display of valid option.
%
% syntax 2: fracGen(mat_p1)
% mat_p1        = user supplied template matrix for the generation. Size 
% should be [n,n,1] or [n,n,n] with integer elements

% additional options: fracGen(...,PERM,NUMp,DO_PLOT)
% PERM          = 0/1: permutate the template matrix at each iteration to 
% create random fractals. Default is 0.
% NUMp          = number of points to be generated. Default is 5000
% DO_PLOT       = 0/1: plot the resulting point distribution and additional 
% details. Default is 1.
%
% % Y.Kamer
% Zurich, 20150428
%
% Coded & Tested on 7.12.0(R2011a)
%%%
%
% Example 1: Generate a regular Merger Sponge fractal
% pts_mat           = fracGen('Menger3D'); 
%
% Example 2: Generate a random 3D monofractal with D=1.58, sampled by 1000 points
% mat_p1            = [0 1; 0 0];
% mat_p1(:,:,2)     = [1 0; 1 0];
% pts_mat           = fracGen(mat_p1,1,1000); 
%
%%%


%% Build in mono/multifractals
% 1D
fracDB.Contor1D             = [0 0 1; 0 0 0; 1 0 0];
fracDB.Contor2D             = [1 0 1; 0 0 0; 1 0 1];
% 2D
fracDB.Vicsek2D             = [0 1 0; 1 1 1; 0 1 0];
fracDB.SierpiTri2D          = [1 0 0 0; 1 1 0 0; 1 0 1 0; 1 1 1 1];
fracDB.SierpiCarp2D         = [1 1 1; 1 0 1; 1 1 1];
fracDB.MultiCarp2D          = [1 1; 1 2];
% 3D
fracDB.OneDin3D(:,:,1)      = [0 1; 0 0];
fracDB.OneDin3D(:,:,2)      = [0 0; 1 0];
fracDB.TwoDin3D(:,:,1)      = [1 1; 0 0];
fracDB.TwoDin3D(:,:,2)      = [0 0; 1 1];
fracDB.Vicsek3D(:,:,1)      = [0 0 0; 0 1 0; 0 0 0];
fracDB.Vicsek3D(:,:,2)      = [0 1 0; 1 1 1; 0 1 0];
fracDB.Vicsek3D(:,:,3)      = [0 0 0; 0 1 0; 0 0 0];
fracDB.Menger3D(:,:,1)      = [1 1 1; 1 0 1; 1 1 1];
fracDB.Menger3D(:,:,2)      = [1 0 1; 0 0 0; 1 0 1];
fracDB.Menger3D(:,:,3)      = [1 1 1; 1 0 1; 1 1 1];
fracDB.MultiCross3D(:,:,1)  = [1 0 2; 0 0 0; 2 0 1];
fracDB.MultiCross3D(:,:,2)  = [0 0 0; 0 1 0; 0 0 0];
fracDB.MultiCross3D(:,:,3)  = [2 0 1; 0 0 0; 1 0 2];
fracDB.Contor3D(:,:,1)      = [1 0 1; 0 0 0; 1 0 1];
fracDB.Contor3D(:,:,2)      = [0 0 0; 0 0 0; 0 0 0];
fracDB.Contor3D(:,:,3)      = [1 0 1; 0 0 0; 1 0 1];

%% Parse input parameters
if(ischar(varargin{1}));
    try
        mat_p1 = fracDB.(varargin{1});
    catch
        fn = fieldnames(fracDB);
        disp('Fractal name should be one of the following:')
        disp(fn);
        pts_out=0;
        return;
    end
else
    mat_p1 = varargin{1};
end

PERM    = 0; 
if(nargin>=2)
    PERM = varargin{2};
end
NUMp    = 5000;
if(nargin>=3)
    NUMp = varargin{3};
end
DO_PLOT = 1;
if(nargin>=4)
    DO_PLOT = varargin{4};
end
sz  = size(mat_p1);
if(sz(1)~=sz(2))
    disp('Template matrix must be a square [n,n,1] or a cube [n,n,n]')
    pts_out=0;
    return;
else
    if(numel(sz)==2)
        MODE_3D=0;
    elseif(numel(sz)==3)
        MODE_3D=1;
        if(sz(3)~=sz(2))
            disp('Template matrix must be a square [n,n,1] or a cube [n,n,n]')
            pts_out=0;
            return;
        end
    else
        disp('Template matrix must be a square [n,n,1] or a cube [n,n,n]')
        pts_out=0;
        return;
    end
end

Mmax    = max(mat_p1(:)); %maximum mass multiplier
Mtot    = sum(mat_p1(:)); %total mass 
iter    = ceil(log(NUMp)/log(Mtot)); %number of iterations to reach NUMp
rndID   = randperm(Mtot^iter); 
rndID   = rndID(1:NUMp);

if(DO_PLOT)
figure;
ax1 = subplot(2, 3, [1]);
ax2 = subplot(2, 3, [2 3 5 6]);
ax3 = subplot(2, 3, [4]);
end
    for i=1:iter
        if(i==1)
            pts_out=fracSingleM([0 0 0 1 0],mat_p1,i,PERM);
        else
            pts_out=fracSingleM(pts_out,mat_p1,i,PERM);
        end
    end
    %% Check for density >1 to fracture further, bypass if input monofractal
    if(Mmax>1)
        IDgt1       = pts_out(:,4)>1;
        nF          = sum(IDgt1); %number of additional frac points
        pts_Frac    = pts_out(IDgt1,:);
        cellOut     = cell(nF,1);
        for f=1:nF
            tmpMass     = pts_Frac(f,4);
            tmpIter     = ceil(log(tmpMass)/log(Mtot/Mmax));
            tmpFracPnt  = [pts_Frac(f,1:3) 1/tmpMass pts_Frac(f,5)];
            for fi=1:tmpIter
                if(fi==1)
                    tmpPts      = fracSingleM(tmpFracPnt,mat_p1,pts_Frac(f,5)+fi,PERM);
                else
                    tmpPts      = fracSingleM(tmpPts,mat_p1,pts_Frac(f,5)+fi,PERM);
                end
            end
            %% Sample the fractal mass without replacement
            P = tmpPts(:,4);
            IDs = nan(tmpMass,1);
            for r = 1:tmpMass
                [~,rndIDD]   = histc(rand,[0 ; cumsum(P(:))] ./ sum(P)) ;
                IDs(r)      = rndIDD;
                P(rndIDD)    = 0;
            end
            tmpPts(:,4)     = 1;
            cellOut{f}      = tmpPts(IDs,:);
        end
        pts_out = [pts_out(~IDgt1,:); cell2mat(cellOut)];
    end
    
    
    if(DO_PLOT)
        %% Plot the generating template matrix
        set(gcf,'CurrentAxes',ax1);
        IDm         = find(mat_p1);
        D           = max(size(mat_p1));
        [Px Py Pz]  = ind2sub(size(mat_p1),IDm);
        
        pts_tmp = [(([Px Py Pz] -0.5 - D/2)/D)];
        if(MODE_3D)
            scatter3(pts_tmp(:,1),pts_tmp(:,2),pts_tmp(:,3),...
                    mat_p1(IDm)*25,mat_p1(IDm),'o','filled',...
                    'markeredgecolor',[1 1 1]*0.5);
        else
            scatter(pts_tmp(:,1),pts_tmp(:,2),...
                    mat_p1(IDm)*25,mat_p1(IDm),'o','filled',...
                    'markeredgecolor',[1 1 1]*0.5);
                
        end
        xlim([-0.5 0.5]);
        ylim([-0.5 0.5]);
        zlim([-0.5 0.5]);
        daspect([1 1 1]);
        if(PERM)
            ttPre = 'Random';
        else
            ttPre = 'Regular';
        end
        if(ischar(varargin{1}))
            title([varargin{1} ' (' ttPre ')'],'fontsize',12);
        else
            title(['User Input (' ttPre ')']','fontsize',12);
        end
        
        set(gcf,'CurrentAxes',ax2);
        scatter3(pts_out(rndID,1),pts_out(rndID,2),pts_out(rndID,3),...
            pts_out(rndID,4)*5,'ok','filled',...
            'markeredgecolor',[1 1 1]*0.2);
        xlim([-0.5 0.5]);
        ylim([-0.5 0.5]);
        zlim([-0.5 0.5]);
        daspect([1 1 1]);
        
        Link = linkprop([ax1, ax2], ...
               {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
        setappdata(gcf, 'StoreTheLink', Link);
        title([ num2str(numel(rndID)) ' points, deepest iteration = ' num2str(max(pts_out(rndID,5))) ],...
            'fontsize',12);
        if(~MODE_3D)
            view(2)
        end
        
        %% Plot the point distribution
        set(gcf,'CurrentAxes',ax3);
        qVec = [-5:0.33:5];
        Dq   = nan(size(qVec));
        for i=1:numel(qVec)
            q       = qVec(i);
            Dq(i)   = (1/(1-q))*(log(sum((mat_p1(mat_p1~=0)/Mtot).^q))/log(size(mat_p1,1)));
        end
        plot(qVec,Dq,'.-k','linewidth',2);
        xlim([-5 5]);
        ylim([floor(min(Dq)) ceil(max(Dq))]);
        pbaspect([1 1 1])
        ylabel('D(q)');
        xlabel('q');
        grid on;
    end
    if(MODE_3D)
        pts_out     = pts_out(rndID,1:3);
    else
        pts_out     = pts_out(rndID,1:2);
    end
end
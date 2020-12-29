%%
clearvars;
close all;
path(pathdef);
%addpath('/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/GHdist/')

addpath('~/gitrepos/SINATRA/Scripts/Data_Generation/GHdist/utils/'); 

path(path, genpath('./utils'));

%%
meshName = 'V13';
raw_mesh_path = '/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/new_aligned_shapesv3/Old_Follivore/';
out_mesh_path = ['/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/new_aligned_shapesv4/' meshName '_2peak_2gpv2bw/'];
touch([out_mesh_path 'mesh/']);

meshNames = getFileNames(raw_mesh_path);
%testIdxes = {[2],[5]};

% testIdxes = {[1,7,3,9,5],[6,2,8,4,10]};
%k = 2

%name = strcat('mesh',int2str(k))
%%
for k = 1:50
    raw_mesh_path = '/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/new_aligned_shapesv3/Old_Follivore/';
    out_mesh_path = ['/Users/brucewang/Documents/new_aligned_shapesv4/' meshName '_5peak_2gp_v' int2str(k) '_bw/'];
    touch([out_mesh_path 'mesh/']);

    meshNames = getFileNames(raw_mesh_path);
    %testIdxes = {[2],[5]};
    idxs1 = randsample(1:10,5, false)
    idxs2 = setdiff(1:10,idxs1)
    testIdxes = {idxs1,idxs2}
    for jj=1:25
        for jjj=1:2
            G = load(['~/gitrepos/SINATRA/GHdist/samples/' meshName '.mat']);
            G = G.G;
            ObsLmkIdx = G.Aux.GTLmks;
            domainMesh = Mesh('VF',G.Aux.UniformizationV,G.F);

            Graw = Mesh('off',[raw_mesh_path 'clean_' meshName '_sas.off']);
            [Area,Center] = Centralize(Graw,'ScaleArea');

            Pts1 = G.V(:, ObsLmkIdx);
            Pts2 = Graw.V(:, ObsLmkIdx);
            [U,~,V] = svd(Pts1*Pts2');
            Rraw = V*U';

            testIdx = ObsLmkIdx(testIdxes{jjj});
            ScalingFactor = ones(G.nV,1);
            ScalingFactor(testIdx) = randi([50,150],size(testIdx,1),size(testIdx,2));
    %         ScalingFactor(testIdx) = randi([50,200],size(testIdx,1),size(testIdx,2));
            A = G.A;
            invD = diag(1./sum(A));
            rwLap = 0.5*invD*A+0.5*eye(size(A,1));

            ScalingSptIdx = zeros(G.nV,1);
            ScalingSptIdx(testIdx) = 1;

            %%% specify a common observer landmark to simulate non-causal
            %%% background fluctuation
    %         ScalingFactor(ObsLmkIdx(commonIdxes)) = randi([100,150],size(commonIdxes,1),size(commonIdxes,2));

            for kk=1:20
                ScalingFactor = rwLap*ScalingFactor;
                ScalingSptIdx = rwLap*ScalingSptIdx;
            end
            ScalingSptIdx = ScalingSptIdx / max(ScalingSptIdx);
    %         ScalingSptIdx = double(ScalingSptIdx > 0);
            csvwrite([out_mesh_path 'gp' num2str(jjj) '_spt.csv'],ScalingSptIdx);

            ScalingFactor = G.F2V*ScalingFactor;
    %         G.ViewFunctionOnMesh(ScalingFactor,struct('mode','native'));
    %         G.ViewFunctionOnMesh(ScalingSptIdx,struct('mode','rb'));

            AffineTransformations = zeros(domainMesh.nF, 3, 2);
            InterpAffineTransformations = zeros(size(AffineTransformations));

            cback = 0;
            for j=1:domainMesh.nF
                localF = domainMesh.F(:,j);
                localV = domainMesh.V(:,localF);
                DomainLocalFrame = [localV(1:2,2)-localV(1:2,1), localV(1:2,3)-localV(1:2,1)];
                localFk = G.F(:,j);
                localVk = G.V(:,localFk);
                kLocalFrame = [localVk(:,2)-localVk(:,1), localVk(:,3)-localVk(:,1)];
                AffineTransformations(j,:,:) = kLocalFrame/DomainLocalFrame;
                InterpAffineTransformations(j,:,:) = ScalingFactor(j) * AffineTransformations(j,:,:);

                for cc=1:cback
                    fprintf('\b');
                end
                cback = fprintf(['%4d/' num2str(domainMesh.nF) ' done.\n'], j);
            end

            [DivX,DivY] = domainMesh.ComputeDivergenceMatrix;
            DivRhoX = DivX*squeeze(InterpAffineTransformations(:,1,1))...
                +DivY*squeeze(InterpAffineTransformations(:,1,2));
            DivRhoY = DivX*squeeze(InterpAffineTransformations(:,2,1))...
                +DivY*squeeze(InterpAffineTransformations(:,2,2));
            DivRhoZ = DivX*squeeze(InterpAffineTransformations(:,3,1))...
                +DivY*squeeze(InterpAffineTransformations(:,3,2));

            %%% reconstruct
            L = domainMesh.ComputeCotanLaplacian;
            reconX = [L;ones(1,domainMesh.nV)]\[DivRhoX;0];
            reconY = [L;ones(1,domainMesh.nV)]\[DivRhoY;0];
            reconZ = [L;ones(1,domainMesh.nV)]\[DivRhoZ;0];

            %% check reconstructed mesh
            reconMesh = Mesh('VF',[reconX';reconY';-reconZ'],domainMesh.F);

            Pts1 = G.V(:, ObsLmkIdx);
            Pts2 = reconMesh.V(:, ObsLmkIdx);
            [U,~,V] = svd(Pts1*Pts2');
            R = V*U';
            reconMesh.V = R*reconMesh.V;
            reconMesh.Centralize('ScaleArea');

            G.V = Rraw*G.V*sqrt(Area)+repmat(Center,1,G.nV);
            reconMesh.V = Rraw*reconMesh.V*sqrt(Area)+repmat(Center,1,reconMesh.nV);

    %         figure;
    %         
    %         h(1) = subplot(1,3,1);
    %         G.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    %         camlight('headlight');
    %         camlight(180,0);
    %         lighting phong;
    %         hold on
    %         scatter3(G.V(1,testIdx),G.V(2,testIdx),G.V(3,testIdx),20,'r','filled');
    %         
    %         h(2) = subplot(1,3,2);
    %         reconMesh.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    %         camlight('headlight');
    %         camlight(180,0);
    %         lighting phong;
    %         hold on
    %         scatter3(reconMesh.V(1,testIdx),reconMesh.V(2,testIdx),reconMesh.V(3,testIdx),20,'r','filled');
    %         
    %         h(3) = subplot(1,3,3);
    %         Graw.draw(struct('FaceColor',[0.9 0.9 0.8],'EdgeColor','none','FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    %         camlight('headlight');
    %         camlight(180,0);
    %         lighting phong;
    %         hold on
    %         scatter3(Graw.V(1,testIdx),Graw.V(2,testIdx),Graw.V(3,testIdx),20,'r','filled');
    %         
    %         Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    %         setappdata(gcf, 'StoreTheLink', Link);
    %         
    %         pause();

            touch([out_mesh_path 'mesh/gp' num2str(jjj) '/']);
            reconMesh.Write([out_mesh_path 'mesh/gp' num2str(jjj) '/' meshName '_' num2str(jj) '.off'], 'off');
        end
    end
end 



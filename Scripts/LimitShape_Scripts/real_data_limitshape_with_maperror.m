clear
close all
clc


addpath('~/gitrepos/SINATRA/Scripts/Data_Generation/GHdist/utils/'); 

addpath('~/gitrepos/LimitShape/Utils/'); 
%% load data: 20 shapes, the first 10 and the second 10 of which consist two classes (interpreted as humans in two different poses)
data_path = '~/Dropbox/SINATRA_Data/cleaned_real_data/';
load('~/Dropbox/SINATRA_Data/cleaned_real_data/Maps_SINATRA_real_data.mat'); %%% this will load the variable "Maps" from the saved file

% set the dimension of the encoded functional maps.
ndim = 60;

% set the error percentage (e.g., 0.25, 0.5, 0.75)
errperc = 0.25;

%%
Grp1 = {'Tarsius'};
Grp2 = {'Microcebus','Mirza'};

Idx1 = [];
Idx2 = [];

for j=1:length(Grp1)
    Idx1 = [Idx1 find(strcmpi(familynames,Grp1{j}))'];
end

for j=1:length(Grp2)
    Idx2 = [Idx2 find(strcmpi(familynames,Grp2{j}))'];
end

Idx_two_groups = [Idx1,Idx2];
Maps = Maps(Idx_two_groups, Idx_two_groups);

Idx1 = 1:length(Idx1);
Idx2 = (length(Idx1)+1):(length(Idx1)+length(Idx2));

%% randomly scramble the maps to generate "artificial error"
for j=1:size(Maps,1)
    for k=1:size(Maps,2)
        if (j==k)
            continue
        else
            %%% scramble the top _errperc_ percent vertices in the map
            n_scram = floor(length(Maps{j,k})*errperc);
            err_seg = Maps{j,k}(1:n_scram);
            Maps{j,k}(1:n_scram) = err_seg(randperm(n_scram));
        end
    end
end

%%
shapes = cell(length(Idx_two_groups), 1);
for j = 1:length(shapes)
    shapes{j} = compute_laplacian_basis(read_off_shape(fullfile(data_path,familynames{Idx_two_groups(j)},filenames{Idx_two_groups(j)})), 60);
    fprintf('%d/%d\n',j,length(shapes));
end

%% Analyze the area-based variability between the two classes. 
% set connectivity graph of the functional maps network: G(i, j) = 1 iff
% Fmaps{i, j} is given. Here we use the dense graph.
G = ones(length(shapes)) - eye(length(shapes));
% convert the p2p maps to functional maps.
Fmaps = convert_fmaps(shapes, Maps, G, ndim);

% compute the consistent latent basis.
U_nn = extract_latent_basis(Fmaps, G); 
% canonicalize the first 40 consistent latent basis.
U = canonical_latent_basis(U_nn, shapes, 40); 

% slect the size of limit shape difference
nCLB = 30; 

% construct limit shape differences.
nshapes_all = length(shapes); 
DA = cell(nshapes_all, 1); 

for i = 1:nshapes_all
    DA{i} = U.area_lsd{i}(1:nCLB, 1:nCLB); 
end

% construct the graphs Gg and Gc. 
% Gg contains edges within the same classes, while Gc contains that across
% different classes.
na = length(Idx1); nb = length(Idx2);
Ag = zeros(nCLB); 
Ac = zeros(nCLB); 
Gg = blkdiag(G(1:na, 1:na), G(na+1:end, na+1:end)); 
Gc = G - Gg; 

for i = 1:nshapes_all
    for j = 1:nshapes_all
        if Gg(i, j) == 1
            Ag = Ag + (DA{i} - DA{j})^2; 
        end
        if Gc(i, j) == 1
            Ac = Ac + (DA{i} - DA{j})^2; 
        end
    end
end

E = Ac/nnz(Gc) - Ag/nnz(Gg); 
[u, v] = eig((E+E')/2); 
[~, id] =sort(diag(v), 'descend'); 
u = u(:, id);             

%% Plot the first few highlighted functions
figure('Name','Area Diff'); 
para.view_angle = [0, 270];
nPerRow = 8;
counter = 1;
for j=1:nPerRow
    fa = (shapes{Idx1(1)}.evecs(:, 1:ndim)*U.cbases{Idx1(1)}(:, 1:nCLB)*u(:, j)).^2;
    fb = (shapes{Idx2(1)}.evecs(:, 1:ndim)*U.cbases{Idx2(1)}(:, 1:nCLB)*u(:, j)).^2;
    h(counter) = subplot(2,nPerRow,j);
    plot_function(shapes{Idx1(1)}, fa, para);
    counter = counter + 1;
    h(counter) = subplot(2,nPerRow,j+nPerRow);
    plot_function(shapes{Idx2(1)}, fb, para);
    counter = counter + 1;
end

cameratoolbar
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);

%%
eigfuncs = cell(1,length(shapes));
mesh_filenames = filenames(Idx_two_groups);
mesh_familynames = familynames(Idx_two_groups);
for j=1:length(shapes)
    eigfuncs{j} = (shapes{j}.evecs(:, 1:ndim)*U.cbases{j}(:, 1:nCLB)*u).^2;
end
save(['./rslt/' num2str(errperc) '_' [Grp1{:}] '_' [Grp2{:}] '_AreaDiff.mat'],...
    'eigfuncs', 'mesh_filenames', 'mesh_familynames', 'Maps');

%% Analyze the extrinsic (i.e., coordinate based) variability across the two classes
DE = cell(nshapes_all, 1); 
% construct extrinsic limit shape difference, please refer to our paper for the
% exact formulation.
for i = 1:nshapes_all
    M = squareform(pdist(shapes{i}.surface.VERT)).^2; 
    M = shapes{i}.A*M*shapes{i}.A; 
    M = diag(sum(M)) - M; 
    
    DE{i} = U.cbases{i}(:, 1:nCLB)'*shapes{i}.evecs'*M*shapes{i}.evecs*U.cbases{i}(:, 1:nCLB); 
end

Ag = zeros(nCLB); 
Ac = zeros(nCLB); 
Gg = blkdiag(G(1:na, 1:na), G(na+1:end, na+1:end)); 
Gc = G - Gg;

for i = 1:nshapes_all
    for j = 1:nshapes_all
        if Gg(i, j) == 1
            Ag = Ag + (DE{i} - DE{j})^2; 
        end
        if Gc(i, j) == 1
            Ac = Ac + (DE{i} - DE{j})^2; 
        end
    end
end

E = Ac/nnz(Gc) - Ag/nnz(Gg); 
[u, v] = eig((E+E')/2); 
[~, id] =sort(diag(v), 'descend'); 
u = u(:, id);             

%% Plot the first few ext highlighted functions
figure('Name','Ext Area Diff'); 
para.view_angle = [0, 270];
nPerRow = 8;
counter = 1;
for j=1:nPerRow
    fa = (shapes{Idx1(1)}.evecs(:, 1:ndim)*U.cbases{Idx1(1)}(:, 1:nCLB)*u(:, j)).^2;
    fb = (shapes{Idx2(1)}.evecs(:, 1:ndim)*U.cbases{Idx2(1)}(:, 1:nCLB)*u(:, j)).^2;
    h(counter) = subplot(2,nPerRow,j);
    plot_function(shapes{Idx1(1)}, fa, para);
    counter = counter + 1;
    h(counter) = subplot(2,nPerRow,j+nPerRow);
    plot_function(shapes{Idx2(1)}, fb, para);
    counter = counter + 1;
end

cameratoolbar
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);

%%
eigfuncs = cell(1,length(shapes));
mesh_filenames = filenames(Idx_two_groups);
mesh_familynames = familynames(Idx_two_groups);
for j=1:length(shapes)
    eigfuncs{j} = (shapes{j}.evecs(:, 1:ndim)*U.cbases{j}(:, 1:nCLB)*u).^2;
end
save(['./rslt/' num2str(errperc) '_' [Grp1{:}] '_' [Grp2{:}] '_ExtAreaDiff.mat'],...
    'eigfuncs', 'mesh_filenames', 'mesh_familynames', 'Maps');


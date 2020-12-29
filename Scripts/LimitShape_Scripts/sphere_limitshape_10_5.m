clear
close all
clc


addpath('~/gitrepos/SINATRA/Scripts/Data_Generation/GHdist/utils/'); 

addpath('~/gitrepos/LimitShape/Utils/'); 
%% load data: 20 shapes, the first 10 and the second 10 of which consist two classes (interpreted as humans in two different poses)
%data_path = './V13_2peak_2gp/';
%errperc = 0.25;
errperc = 0;
%nvertex = 5131;
nvertex = 642;

nfuncs = 15;
times = [];
for a = 1:100
    tic

    data_path = ['/Users/brucewang/Documents/spheres10_5_2020-05-15/v' num2str(a) '/']; 

   % NumPerGroup = 25;
    NumPerGroup = 50;

    NumGroups = 2;

    shapes = cell(NumPerGroup*NumGroups, 1);
    for i = 1:NumPerGroup
        %shapes{i} = compute_laplacian_basis(read_off_shape([data_path 'mesh/gp1/V13_' num2str(i) '.off']), 60); 
        %shapes{i+NumPerGroup} = compute_laplacian_basis(read_off_shape([data_path 'mesh/gp2/V13_' num2str(i) '.off']), 60); 
        shapes{i} = compute_laplacian_basis(read_off_shape([data_path 'v1/sphere_' num2str(i) '.off']), 60); 
        shapes{i+NumPerGroup} = compute_laplacian_basis(read_off_shape([data_path 'v2/sphere_' num2str(i) '.off']), 60); 
      
      %For Caricatures
  %    shapes{i} = compute_laplacian_basis(read_off_shape([data_path 'gp1/V13_' num2str(i) '.off']), 60); 
    %  shapes{i+NumPerGroup} = compute_laplacian_basis(read_off_shape([data_path 'gp2/V13_' num2str(i) '.off']), 60); 

    end


    % Analyze the area-based variability between the two classes. 

    % set the dimension of the encoded functional maps.
    ndim = 60;

    % set maps across the shapes of interest, here all the shapes are in identity
    % to each other. 
    Maps = cell(NumPerGroup*NumGroups, 1); 
    for i = 1:NumPerGroup*NumGroups
        for j = 1:NumPerGroup*NumGroups
            Maps{i, j} = 1:nvertex;
        end
    end
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

    % set connectivity graph of the functional maps network: G(i, j) = 1 iff
    % Fmaps{i, j} is given. Here we use the dense graph.
    G = ones(NumPerGroup*NumGroups) - eye(NumPerGroup*NumGroups);
    % convert the p2p maps to functional maps.
    Fmaps = convert_fmaps(shapes, Maps, G, ndim);

    % compute the consistent latent basis.
    U_nn = extract_latent_basis(Fmaps, G); 
    % canonicalize the first 40 consistent latent basis.
    U = canonical_latent_basis(U_nn, shapes, 40); 

    % slect the size of limit shape difference
    nCLB = nfuncs; 

    % construct limit shape differences.
    nshapes_all = length(shapes); 
    DA = cell(nshapes_all, 1); 

    for i = 1:nshapes_all
        DA{i} = U.area_lsd{i}(1:nCLB, 1:nCLB); 
    end

    % construct the graphs Gg and Gc. 
    % Gg contains edges within the same classes, while Gc contains that across
    % different classes.
    na = NumPerGroup; nb = NumPerGroup;
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
    
   % fb = zeros(642,1)
 %   fa = zeros(642,1)
    
    for index = 1:NumPerGroup
        fa_table = zeros(nvertex,nfuncs);
        fb_table = zeros(nvertex,nfuncs);
        for i = 1:nfuncs
          %  fb = fb + ((shapes{NumPerGroup+index}.evecs(:, 1:ndim)*U.cbases{NumPerGroup+index}(:, 1:nCLB)*u(:, i)).^2)/i
          %  fa = fa + ((shapes{index}.evecs(:, 1:ndim)*U.cbases{index}(:, 1:nCLB)*u(:, i)).^2)/i
            fa_table(:,i) = ((shapes{index}.evecs(:, 1:ndim)*U.cbases{index}(:, 1:nCLB)*u(:, i)).^2);
            fb_table(:,i) = ((shapes{NumPerGroup+index}.evecs(:, 1:ndim)*U.cbases{NumPerGroup+index}(:, 1:nCLB)*u(:, i)).^2);
        end 
      %  csvwrite([data_path 'gp1/map_' num2str(index) '.csv'],fa_table);
       % csvwrite([data_path 'gp2/map_' num2str(index) '.csv'],fb_table);
      %  csvwrite([data_path 'v1/map_' num2str(index) '.csv'],fa_table);
      %  csvwrite([data_path 'v2/map_' num2str(index) '.csv'],fb_table);
    end 
    timeElapsed = toc;
    times = [times, timeElapsed]

end 
csvwrite(['~/gitrepos/SINATRA/Scripts/Data/' num2str(nfuncs) 'spheres.csv'], times)
%% Plot the first highlighted function
figure; 
para.view_angle = [0, 270]; 

fa = (shapes{1}.evecs(:, 1:ndim)*U.cbases{1}(:, 1:nCLB)*u(:, 1)).^2;
fb = (shapes{NumPerGroup+1}.evecs(:, 1:ndim)*U.cbases{NumPerGroup+1}(:, 1:nCLB)*u(:, 1)).^2; 
h(1) = subplot(121); plot_function(shapes{1}, fa, para); title('Area Diff First: Class A'); 
h(2) = subplot(122); plot_function(shapes{NumPerGroup+1}, fb, para); title('Area Diff First: Class B');
cameratoolbar
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
 
u = u(:, id);             

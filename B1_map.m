%% Clear workspace
clearvars; clc; close all;

%% Select folder with NIfTI files
addpath('/Users/cbd/Downloads/Matfiles_New/ColorMaps');
disp('Select B1 map folder');
DIR = uigetdir(pwd, 'Select B1 map containing NIfTI files (.nii)');
niiFiles = dir(fullfile(DIR, '*.nii'));

if isempty(niiFiles)
    error('No NIfTI files found in this folder!');
end

%% Load first NIfTI file
filePath = fullfile(DIR, niiFiles(1).name);
niiVolume = double(niftiread(filePath));  % convert to double for calculations
info = niftiinfo(filePath);

% Check dimensions
dims = size(niiVolume);  % e.g., [128, 256, 1, 4] or [128,256,4]
xDim = dims(1); 
yDim = dims(2); 
zDim = dims(3); 
if numel(dims) >= 3
    zDim = dims(3);
else
    zDim = 1;
end
if numel(dims) == 4
    nVol = dims(4);
else
    nVol = 1;
end

disp(['Loaded NIfTI with dimensions: ', num2str(dims)]);
if nVol > 1
    v = 2;  % GE B1 maps often stored in volume 2
    S = niiVolume(:,:,1,v);
else
    v = 1;
    S = niiVolume(:,:,1);
end

%% Automatically detect the phantom outline instead of drawing manually
S = niiVolume(:,:,:,2); 

% Step 1: Smooth the image to reduce noise
Sm = imgaussfilt(S, 4);

% Step 2: Normalize image to [0,1]
Smn = mat2gray(Sm);

% Step 3: Threshold (Otsu)
bw_phantom = imbinarize(Smn);

% Step 4: Keep only the largest connected component (phantom box)
CCp = bwconncomp(bw_phantom);
numPix = cellfun(@numel, CCp.PixelIdxList);
[~, idxMax] = max(numPix);
phantom_outline = false(size(bw_phantom));
phantom_outline(CCp.PixelIdxList{idxMax}) = true;

% Step 5: Fill any holes inside phantom
phantom_outline = imfill(phantom_outline, "holes");

% Optional: Smooth outline
phantom_outline = imopen(phantom_outline, strel('disk', 5));

% Show detected phantom outline
figure; 
imshowpair(S, phantom_outline, 'montage');
title('Auto-detected Phantom Outline (right)');

% Overlay outline in red on original
figure; 
imshow(S, [], 'InitialMagnification', 'fit'); colormap("gray"); hold on;
visboundaries(phantom_outline, 'Color','r','LineWidth',2); 
title('Automatic Phantom Outline');
hold off;
%% Ask user if the automatic outline is acceptable
choice = questdlg('Is the automatic phantom outline acceptable?', ...
    'Outline Verification', ...
    'Yes, continue','No, manually adjust','No, manually adjust');

% If manual adjustment is needed
if strcmp(choice, 'No, manually adjust')
    fprintf('Manual outline adjustment selected.\n');
    
    % Create figure for manual drawing
    fig_manual = figure;
    imshow(S, [], 'InitialMagnification', 'fit'); 
    colormap("gray");
    title('Draw phantom outline manually. Double-click to finish.');
    
    % Allow user to draw polygon
    h = drawpolygon('Color','r','LineWidth',2);
    
    % Wait for user to double-click
    wait(h);
    
    % Create mask from the drawn polygon
    phantom_outline = createMask(h);
    
    % Close the manual drawing figure
    close(fig_manual);
    
    % Show the final manual outline
    figure;
    imshow(S, [], 'InitialMagnification', 'fit'); colormap("gray"); hold on;
    visboundaries(phantom_outline, 'Color','r','LineWidth',2);
    title('Final Manual Phantom Outline');
    hold off;
    
else
    fprintf('Automatic outline accepted. Continuing...\n');
end

%% --- Normalize B1 map to phantom mean ---
B1map = squeeze(niiVolume(:,:,1,v));      % original B1 map (scanner units)
B1_phantom = B1map(phantom_outline);      % values inside phantom

% Normalize: mean inside phantom = 100%
B1map_percent = (B1map ./ mean(B1_phantom)) * 100;

% Mask outside phantom
B1map_percent(~phantom_outline) = NaN;

%% --- Compute statistics inside phantom ---
mean_B1_percent = mean(B1map_percent(phantom_outline),'omitnan');
std_B1_percent  = std(B1map_percent(phantom_outline),'omitnan');

%% --- Plot first slice as example ---
load('T1cm.mat')
load('T2cm.mat')
load('differenceMaps.mat')


figure;
ax1 = axes;
imagesc(niiVolume(:,:,1,2))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, B1map_percent, 'alphadata', phantom_outline);
axis image off;
colormap(ax2,T2colormap);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('B1 map');
colorbar;
set(gca, 'FontSize', 24);
saveas(gcf, 'B1map_0.tiff');
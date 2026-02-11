%% Clear workspace
clearvars; clc; close all;

%% Select folder containing DICOM files
disp('Select folder containing DICOM files (B0 map)');
DIR = uigetdir(pwd, 'Select folder with DICOM files');
if DIR == 0
    error('No folder selected!');
end

%% Get list of DICOM files
dicomFiles = dir(fullfile(DIR, '*.dcm'));
if isempty(dicomFiles)
    error('No DICOM files found in this folder!');
end

%% Read DICOM images
numFiles = numel(dicomFiles);
info = dicominfo(fullfile(DIR, dicomFiles(1).name));
dims = [info.Rows, info.Columns, numFiles];
B0_Hz = zeros(dims);

for k = 1:numFiles
    fname = fullfile(DIR, dicomFiles(k).name);
    B0_Hz(:,:,k) = double(dicomread(fname));
end

disp(['Loaded DICOM volume with dimensions: ', num2str(size(B0_Hz))]);

%% If only one slice, take the first slice
B0_Hz = B0_Hz(:,:,1);

%% Automatically detect the phantom outline
Sm = imgaussfilt(B0_Hz, 2);          % Smooth
Smn = mat2gray(Sm);                  % Normalize
bw_phantom = imbinarize(Smn);        % Threshold (Otsu)
CCp = bwconncomp(bw_phantom);        % Connected components
numPix = cellfun(@numel, CCp.PixelIdxList);
[~, idxMax] = max(numPix);
phantom_outline = false(size(bw_phantom));
phantom_outline(CCp.PixelIdxList{idxMax}) = true;
phantom_outline = imfill(phantom_outline,'holes');       % Fill holes
phantom_outline = imopen(phantom_outline, strel('disk',3)); % Optional smooth

%% Show overlay
figure; imshow(B0_Hz, [], 'InitialMagnification','fit'); hold on;
visboundaries(phantom_outline,'Color','r','LineWidth',2);
title('Phantom Outline on B0 Map'); hold off;

%% Convert from Hz to ppm
B0_field = 3;          % Tesla
f0_MHz = 42.577*B0_field;  % Larmor frequency
B0_ppm = B0_Hz ./ f0_MHz;

%% Clip to [-1, +1] ppm
B0_ppm(B0_ppm < -1) = -1;
B0_ppm(B0_ppm >  1) = 1;

%% Display B0 ppm map
figure;
imagesc(B0_ppm, [-1 1]);
axis image off;
colormap(jet); colorbar;
title(sprintf('B0 Map (ppm, %.1fT, range -1 to +1)', B0_field));
%% --- Plot first slice as example ---
% load('T1cm.mat')
% load('T2cm.mat')
% load('differenceMaps.mat')
% figure;
% ax1 = axes;
% imagesc(niiVolume(:,:,1,1))
% axis image off;
% colormap(ax1, 'gray');
% ax2 = axes;
% imagesc(ax2, B1map_percent, 'alphadata', phantom_outline);
% axis image off;
% colormap(ax2,cm);
% ax2.Visible = 'off';
% linkprop([ax1 ax2],'Position');
% title('B1 map');
% colorbar;
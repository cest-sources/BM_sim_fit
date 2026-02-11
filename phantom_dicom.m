clear all; clc; close all;
folderPath = '/Users/cbd/Downloads/CESTPhantomCode/flywheel/RMRCBD061025/20250610_3549/15 - SSFSE, CW, DL';
addpath('/Users/cbd/Downloads/Matfiles_New/ColorMaps');
load('T2cm.mat');
load('T1cm.mat');
load('differenceMaps.mat');
% Get list of DICOM files
files = dir(fullfile(folderPath, '*.dcm'));

% Sort files by name (important if slices are ordered in filenames)
[~, idx] = sort({files.name});
files = files(idx);

% Preallocate cell array for images
numFiles = length(files);
imgs = cell(1, numFiles);

% Read all DICOM images
for k = 1:numFiles
    fname = fullfile(folderPath, files(k).name);
    imgs{k} = dicomread(fname);
end

% Convert to 4D array for montage
imgs4D = cat(4, imgs{:});

% Display all 63 slices in one figure
figure;
montage(imgs4D, 'DisplayRange', []);
colormap(T1colormap);colorbar;

% figure;
% montage(imgs4D, 'DisplayRange', []);
% colormap(T1colormap);
% 
% figure;
% montage(imgs4D, 'DisplayRange', []);
% colormap(cm);
% 
% figure;
% montage(imgs4D, 'DisplayRange', []);
% colormap(cm1);
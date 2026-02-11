% calcMTRmap: Calculates the MTR asymmetry map using the z-spectroscopic 
% imaging data. If a specific ppm value is specified, only the MTR
% asymmetry value pertaining to that value will be returned; otherwise, all
% MTR asymmetry values for each positive-negative offset pair will be
% returned. MTR asymmetry is Z_negativePPM - Z_positivePPM.
%
%   INPUTS:
%       zspecImgData    -   Array containing z-spectroscopic (imaging) data.
%                           Last dimension in the array is the
%                           spectroscopic dimension.
%       ppmOffsets      -   Numeric vector containing the saturation 
%                           offsets (in ppm) matching the spectroscopic
%                           dimension of zspecImgData
%       sel_ppm         -   (optional) Number pertaining to ppm value to 
%                           use to calculate MTR asymmetry map; otherwise,
%                           MTR asymmetry is calculated across all ppm
%                           values
%
%   OUTPUTS:
%       MTRmap          -   Array containing MTR asymmetry image or spectra. 
%                           All but the last dimension are identical to 
%                           zspecImgData. The last dimension will be 1 if
%                           sel_ppm was specified; otherwise it will be the
%                           length of positive-negative offset pairs
%                           identified
%       MTRppm          -   Vector of positive ppm values pertaining to the
%                           (now-reordered) MTRmap spectral dimension
%       sel_ppm_true    -   (if sel_ppm specified as input) Number
%                           pertaining to actual ppm value used to 
%                           calculate MTR asymmetry map. If sel_ppm not
%                           found in ppmOffsets, it will return the 1st 
%                           value larger than sel_ppm in ppmOffsets
%
function [MTRmap,MTRppm,sel_ppm_true]=calcMTRmap(zspecImgData,ppmOffsets,sel_ppm)

% Detect whether input sel_ppm was specified
if nargin<3
    selppmflg=false;
else
    selppmflg=true;
end

zsIDsize=size(zspecImgData);
if length(zsIDsize)>2 
    % Reshape zspecImgData to have only 2 dimensions
    zspecImgData=reshape(zspecImgData,[],zsIDsize(end));
elseif length(zsIDsize)==1
    % Handle 1D case (single spectrum)
    zspecImgData=reshape(zspecImgData,1,[]);
    zsIDsize=[1 length(zspecImgData)];
end

% Reorder both ppmOffsets and zspecImgData (in spectroscopic dimension) to
% be increasing in offset
[ppmOffsets,sortIdx]=sort(ppmOffsets,'ascend');
zspecImgData=zspecImgData(:,sortIdx);

% Determine the z-spec values for the positive and negative ppm values
if selppmflg
    % Find the indices pertaining to positive and negative sel_ppm
    selPpmPosIdx=find(abs(sel_ppm-ppmOffsets)==min(abs(sel_ppm-ppmOffsets)));
    if length(selPpmPosIdx)>1 %select the larger abs value
        selPpmPosIdx=selPpmPosIdx(end);
    end
    selPpmNegIdx=find(abs(-sel_ppm-ppmOffsets)==min(abs(-sel_ppm-ppmOffsets)));
    if length(selPpmNegIdx)>1 %select the larger abs value
        selPpmNegIdx=selPpmNegIdx(1);
    end

    zImgPos=zspecImgData(:,selPpmPosIdx);
    zImgNeg=zspecImgData(:,selPpmNegIdx);

    sel_ppm_true=ppmOffsets(selPpmPosIdx);
else
    % Go through the positive offsets and find the negative offset that is
    % closest in magnitude to each
    ppmPosIdx=find(ppmOffsets>0);
    ppmPosIdxPaired=ppmPosIdx;
    ppmNegIdxPaired=zeros(size(ppmPosIdxPaired));
    
    % Get all negative ppm values and their indices in the full array
    ppmNegIdx=find(ppmOffsets<0);
    ppmNegVals=ppmOffsets(ppmNegIdx);
    
    for ii=1:length(ppmPosIdxPaired)
        % Find the negative offset closest in magnitude to current positive offset
        targetVal = -ppmOffsets(ppmPosIdxPaired(ii));
        [~, minIdx] = min(abs(ppmNegVals - targetVal));
        
        % Store the index in the original ppmOffsets array
        ppmNegIdxPaired(ii) = ppmNegIdx(minIdx);
    end
    
    % Finally, ID the positive and negative offset Z-values
    zImgPos=zspecImgData(:,ppmPosIdxPaired);
    zImgNeg=zspecImgData(:,ppmNegIdxPaired);

    sel_ppm_true=[]; %no value returned if no sel_ppm specified
end

% Calculate the MTR asymmetry (and corresponding positive ppm offsets)
MTRmap=squeeze(zImgNeg-zImgPos);
MTRppm=ppmOffsets(ppmOffsets>0);

% Reshape MTRmap to match input zspecImgData (except the last dimension!)
MTRmsize=size(MTRmap);
if length(zsIDsize)>2    
    MTRmap=reshape(MTRmap,[zsIDsize(1:end-1) MTRmsize(end)]);
end
end
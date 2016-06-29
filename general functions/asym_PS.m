function [x_asym, ASYM] = asym_PS(x_zspec,Z)
%Calculates asymmetrie of z-spectra
%
%   input:  x_zspec (vector)
%           Z (1D z-spectrum or 4D or 3D zspec stack)
%
%   output: x_asym (vector)
%           ASYM (asym vector or 4D asym stack )

% get dimension of Z stack
Z_dim = numel(size(Z));

% % reshape old 3D-stacks into 4D stack [(x,y,offset) -> (x,y,slice,offset)]
if Z_dim == 3
        Z = reshape(Z,[size(Z,1),size(Z,2),1,size(Z,3)]);   
end

% new x vector
int1    = x_zspec(2:end);
int2    = x_zspec(1:end-1);
step    = abs(min(int1-int2));
offset  = max([abs(min(x_zspec)) abs(max(x_zspec))]);
x_zspec_int = -offset:step:offset;
int3=fix(numel(x_zspec_int)/2):-1:1;
int4=fix(numel(x_zspec_int)/2)+2:numel(x_zspec_int);
x_asym  = 0:step:offset;

if Z_dim == 2
    
    % read z-spectrum
    zspec = Z;
    
    % interpolate zspec data
    y_zspec_int = spline(x_zspec,zspec,x_zspec_int);

    % calculate y_asym vector
    ASYM(1) = 0;
    ASYM(2:numel(int4)+1)=(y_zspec_int(int3)-y_zspec_int(int4));

else

    % x,y,slice loop
    for ii = 1:size(Z,1)
        for jj = 1:size(Z,2)
            for kk = 1:size(Z,3)

                % read z-spectrum
                zspec = squeeze(Z(ii,jj,kk,:));
                
                if (sum(isnan(zspec)) == 0) 

                    % interpolate zspec data
                    y_zspec_int = spline(x_zspec,zspec,x_zspec_int);

                    % calculate y_asym vector
                    y_asym(1) = 0;
                    y_asym(2:numel(int4)+1)=(y_zspec_int(int3)-y_zspec_int(int4));

                    % write into ASYM stack
                    ASYM(ii,jj,kk,:) = y_asym;
                else
                    ASYM(ii,jj,kk,:) = NaN(1,numel(int4)+1);
                end

            end
        end
    end
end

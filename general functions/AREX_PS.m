function [x_arex, AREX] = AREX_PS(x_zspec,Z,R1A)
%Calculates AREX of z-spectra
%
%   input:  Z (1D z-spectrum or 3D/4D zspec stack)
%           R1A (R1A-value/2D-matrix/3D-matrix)
%           
%   output: AREX (vector or 4D stack)

if nargin < 3
    R1A = 1;
    disp('No R1A-value given - calculating MTR_Rex instead of AREX');
end

% get dimension of Z stack
Z_dim = numel(size(Z));
R1_dim = numel(size(R1A));

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
x_arex  = 0:step:offset;

if Z_dim == 2 
    if size(R1A,1) == 1 && size(R1A,2) == 1  
        % read z-spectrum
        zspec = Z;
        % interpolate zspec data
        y_zspec_int = spline(x_zspec,zspec,x_zspec_int);
        % calculate AREX-vector
        AREX(1)=0;
        AREX(2:numel(int4)+1)=(1./y_zspec_int(int4)-1./y_zspec_int(int3))*R1A;
    else
        error('Dimensions of Z and R1A does not match')
    end
else
    % x,y,slice loop
    for ii = 1:size(Z,1)
        for jj = 1:size(Z,2)
            for kk = 1:size(Z,3)

                % read z-spectrum
                zspec = squeeze(Z(ii,jj,kk,:));
                % interpolate zspec data
                y_zspec_int = spline(x_zspec,zspec,x_zspec_int);
                % read R1A
                if R1_dim == 2
                    if size(R1A,1) == 1 && size(R1A,2) == 1
                        r1a = R1A;
                    else
                        r1a = R1A(ii,jj);
                    end
                elseif R1_dim == 3
                    r1a = R1A(ii,jj,kk);
                else
                    error('Could not read R1A')
                end
                
                if (sum(isnan(zspec)) == 0)
                    % calculate AREX-vector
                    y_arex(1)=0;
                    y_arex(2:numel(int4)+1)=(1./y_zspec_int(int4)-1./y_zspec_int(int3))*r1a;

                    % write into AREX stack
                    AREX(ii,jj,kk,:) = y_arex;
                else
                    AREX(ii,jj,kk,:) = NaN(1,numel(x_zspec_int));
                end
            end
        end
    end  
end

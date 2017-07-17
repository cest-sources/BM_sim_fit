function Zrow = b0correct_run(Zrow,runstr)
%  Ztab('ucl1') = norm_run(Ztab('ucl1'),1,0)
    
    w=Zrow.offset{1};
    wfit = interp1(1:numel(w),w,linspace(1,numel(w),numel(w)*20))'; % interpolated w-axes 
    subtab=Zrow{:,runstr};
    
    % normalize
    II=1:numel(subtab); % for all run entries
    
    for ii=II
        if ~isempty(subtab{ii})
            Z=subtab{ii};
            
            Z_splined = interp1(w,bfilt(Z)',wfit,'spline');
           
            [minval,minind]=min(Z_splined);
            dB0=wfit(minind,1);

            if abs(dB0)>1 %frage ab: aktuelle B0 verschiebung
               warning('DANGER in B0 correction - dB0 larger than 1 ppm') 
            end;
            
            Z_lint = interp1(w,Z',wfit,'linear');
            www=wfit-wfit(minind);
            subtab{ii}= interp1(www,Z_lint,w,'linear','extrap');

        end;
    end;
    
    Zrow{:,runstr}=subtab;








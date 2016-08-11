
%Ztab batch

%% load from a array
%% loadQUESP
C=mat2cell(squeeze(I89(:,:,1:11)),[77],[1 1 1 1 1 1 1 1 1 1 1])

Ztab.B1_run('ucl3',1:11)=C

%% loadQUEST
C=mat2cell(squeeze(I89(:,:,12:18)),[77],[1 1 1 1 1 1 1]);

Ztab.tsat_run('ucl2p',1:8) = {7.5, 6 , 5 , 4.0 , 3.0 , 2.0 , 1.0, 0};
Ztab.tsat_run('ucl2',2:8) = C

%% normalize B1run
Ztab('ucl3','B1_run') = norm_run(Ztab('ucl3','B1_run'),1,0)

%% normalize tsatrun
Ztab('ucl2','tsat_run') = norm_run(Ztab('ucl2','tsat_run'),5,0)
%%

create_Ztab();


%%
[xxx, xZspec, Z, varyval, vary, Ptab]=plot_tab(Ztab,'ucl1','B1_run')
%%
[xxx, xZspec, Z, varyval, vary, Ptab]=plot_tab(Ztab,'ucl2','tsat_run')

%% import from zrun
for ii=1:7
Ztab.B1_run{3,ii}=1
Ztab.B1_run{4,ii}=ZRUN{1, ii}{1, 1}.Zmean;
end;

Ztab.offset{4}=ZRUN{1, 1}{1, 1}.x;
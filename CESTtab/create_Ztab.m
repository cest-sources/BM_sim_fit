
offset=-5:5;
Z= [1 0.9 0.8 0.7 0.5 0 0.5 0.7 0.8 0.9 1];

P.readout='GRE'; %readout type
P.readout_flipangle=10; % readout excitation flipangle

P.shape='SPINLOCK';
P.B1=1;  % [µT] standard saturation pulse amplitude (block pulse!)
P.tsat=6;  % [2] standard saturation pulse duration
P.Trec=6; % [s] recover time after readout, before saturation (most probably TR-tsat)
P.B0=7;  % [T] scanner static field strength


%% create table
Ztab=table('RowNames',{'x','ucl1','ucl1','exp3'});
Ztab.Var1{1}=[];
Ztab.Var2{1} =[];
Ztab.Var3{1}=[];
Ztab.Var4{1,1}=1.5
Ztab.Var4{1,2}=1.5;
Ztab.Var4{1,3}=2.25;
Ztab.Var4{1,4}=2.25;
Ztab.Var5{1,1}=1.5;
Ztab.Var5{1,2}=1.5;
Ztab.Var5{1,3}=2.25;
Ztab.Var5{1,4}=2.25;
Ztab.Properties.VariableNames{1} = 'exp';
Ztab.Properties.VariableNames{2} = 'P';
Ztab.Properties.VariableNames{3} = 'offset';
Ztab.Properties.VariableNames{4} = 'B1_run';
Ztab.Properties.VariableNames{5} = 'tsat_run';

%% fill table for QUESP
% used B1 values

Ztab.B1_run{1,1} = 0.5;
Ztab.B1_run{1,2} = 0.75;
Ztab.B1_run{1,3} = 1;
Ztab.B1_run{1,4} = 1.5;

expNR=1;
Ztab.exp{expNR+1} = 'pH=7, buffer=5mM';
Ztab.P{expNR+1} = P;
Ztab.offset{expNR+1,1} = offset;
Ztab.B1_run{expNR+1,1} = Z;
Ztab.B1_run{expNR+1,2} = 0.9*Z;
Ztab.B1_run{expNR+1,3} = 0.8*Z;
Ztab.B1_run{expNR+1,4} = 0.7*Z;
% plot
for ii=1:4
figure(42), plot(Ztab.offset{2},Ztab.B1_run{2,ii}); hold on;
end
title('QUESP data');
Ztab
Ztab{:,{'B1_run'}}

%% fill table for QUEST
% used B1 values
Ztab.tsat_run{1,1} = 0.5;
Ztab.tsat_run{1,2} = 0.75;
Ztab.tsat_run{1,3} = 1;
Ztab.tsat_run{1,4} = 1.5;

Ztab.offset{2,1} = offset;
Ztab.tsat_run{2,1} = Z;
Ztab.tsat_run{2,2} = 0.9*Z;
Ztab.tsat_run{2,3} = 0.8*Z;
Ztab.tsat_run{2,4} = 0.7*Z;
% plot
for ii=1:4
figure(43), plot(Ztab.offset{2},Ztab.tsat_run{2,ii}); hold on; 
end
title('QUEST data');

Ztab{:,{'tsat_run'}}


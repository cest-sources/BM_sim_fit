function  [w_x, Z_x, P, rowname, Ztab]=LOAD_xls_qCEST_2_Ztab(sheets)
%%
[file, path] =uigetfile('*');
if nargin<1
[~, sheets]=xlsfinfo([path file]);
end;
% sheets = {    '50mM PB'    '20mM PB'    '10mM PB'    '5mM PB' };
% Ztab % warning is normal
Ztab=table('RowNames',sheets);
Ztab.Var1{1}=[];
Ztab.Var2{1}=[];
Ztab.Var3{1}=[];
Ztab.Var4{1,6}=[];

Ztab.Properties.VariableNames{1} = 'exp';
Ztab.Properties.VariableNames{2} = 'P';
Ztab.Properties.VariableNames{3} = 'offset';
Ztab.Properties.VariableNames{4} = 'B1_run';
clc
%


for JJ=1:numel(sheets)
tablestr=sheets{JJ};
%%%% exp name
[~, ~, raw] = xlsread([path file],tablestr,'B1');
expname=raw{1};

Ztab.Properties.RowNames{JJ} = [expname(1:4) sheets{JJ}] ;
rowname=Ztab.Properties.RowNames{JJ};
%%%% exp params
[~, ~, raw] = xlsread([path file],tablestr,'B2:B15');

params = reshape([raw{:}],size(raw));

P.FREQ = params(1); % e.g. 7 * 42.5764; [MHz] scanner static field strength
P.Trec=params(2); 
P.tp=params(3); 
P.R1A=1/params(4);
P.R2A=1/params(5);
P.cmM=params(6);
P.M0=params(7);
P.Zi=params(8);
P.normalized=params(9);
if params(9)==0 P.normalized=[]; end;
NEX=params(10);
NB1=params(11);
P.pulsed=params(12);
P.n=params(13); % number of stauration pulses
P.DC=params(14); %duty cycle
P.B1=1; % [µT] standard saturation pulse amplitude (block pulse!) (need to add the train of gaussian pulses)

[~, ~, raw] = xlsread([path file],tablestr,'B16'); % pulse shape
P.shape= raw{1};

%%%% load Z

loadstr=sprintf('B24:%s%d',char(double('A')+NB1),23+NEX); % this collects the right data
[~, ~, raw] = xlsread([path file],tablestr,loadstr);
Z_x = reshape([raw{:}],size(raw));
Z_xxx = reshape(Z_x,size(Z_x,1)*size(Z_x,2),1);

Z_xxx=Z_xxx./P.M0;
Z_x=Z_x./P.M0;

%%%% load offsets
[~, ~, raw] = xlsread([path file],tablestr,sprintf('A24:A%d',23+NEX));
P.xZspec = reshape([raw{:}],size(raw));
w_x=P.xZspec;

%%%% load B1
loadstr=sprintf('B23:%s23',char(double('A')+NB1)); % this collects the right data

[~, ~, raw] = xlsread([path file],tablestr,loadstr);
P.varyval=[];
for jjj=1:numel(raw);
    P.varyval =[ P.varyval raw{jjj}];
end;
P.vary{1}='B1';

% Clear temporary variables
% cleavars raw BlocheqGoran;



row=JJ;
for ii=1:numel(P.varyval)
Ztab.B1_run(row,ii) = {Z_x(:,ii)};
end;
Ztab.exp{row} = expname;
Ztab.P{row} = P;
Ztab.offset{row}=P.xZspec;

% normalization and exclusion
% Ztab(row,:) = norm_run(Ztab(row,:),'B1_run',P.normalized) ;
 
figure,
plot_tab(Ztab,Ztab.Properties.RowNames{JJ},'B1_run');

end;

%%%% extend B1run to match others
Ztab.B1_run{1,11} = []




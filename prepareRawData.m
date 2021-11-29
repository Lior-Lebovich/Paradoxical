%% Read raw data and save to csv table:

% raw data is available in all_mice_simplified_data_tables.mat - matlab 
% data file. 
% mouse_data_tables is a 1x6 cell array, where mouse_data_tables{m} is a
% the data table of mouse m. Each row corresponds to an individual trial.
% Variables in each mouse table:
% res: response string: 'r' = right lick; 'l' = left lick; 'NA' = invalid
% trial (see methods).
% targetCondition: 'A' = target A; 'B' = target B; 'N' = no target.
% bgN: number of background odors (0-5 when target A or B are on; 1-6 when 
% no target is on).
% RT: reaction-time in ms.
% normRT: normalized reaction-time in ms (see methods).

% read raw data:

load('all_mice_simplified_data_tables.mat');

tarTypes = {'A','B','N'}; 
minTarN = [1,1,0];

% run over mice:
for m = 1:6
    % read mouse table:
    mouseTableRaw = mouse_data_tables{m};
    % omit invalid trials:
    mouseTable = mouseTableRaw( ~strcmp( mouseTableRaw.res, 'NA' ) , : );
    % add mouse id (0-5):
    mouseTable.subj_idx = (m-1) * ones( size(mouseTable,1), 1 );
    % define 18 conditions (6 possible #BG-odors X 3 possible target types):
    mouseTable.stim = mouseTable.targetCondition + ...
        num2str( mouseTable.bgN );
    % save RT in sec:
    mouseTable.rtNorm = mouseTable.normRT / 1000;
    % save logical response (1 = "target present" choice; 0 = "target absent" choice)
    mouseTable.response = strcmp( mouseTable.res, 'r' );    
    intermData.(['m' num2str(m)]) = mouseTable;
    clear mouseTable;
end

allMiceTable = intermData.m1( 1:0, : );
for m = 1:6
    allMiceTable = [allMiceTable; intermData.(['m' num2str(m)])];
end

% order columns and omit irrelevant variables:
miceTableCond_ordCol_RTNorm = removevars( allMiceTable, ...
    {'res','targetCondition','bgN','RT','normRT'} );

% save tables as CSV - normalized RTs:
miceTableCond_ordCol_RTNorm.Properties.VariableNames{'rtNorm'} = 'rt';
writetable( miceTableCond_ordCol_RTNorm, 'goodAllMice_id_stim_rtNorm_choiceQ_0.csv' );



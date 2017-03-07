function [data] = readincortex(indirRL,indirLR,hemisphere)
% Benjamin Risk
% This function extracts the specified hemisphere from the cifti file
% and converts to a matlab matrix. The data from session 1 (RL phase 
% encoding) and session 2 (LR phase encoding) are concatenated.

hemisphere = validatestring(hemisphere, {'CORTEX_RIGHT','CORTEX_LEFT'});
side = 'R';
if strcmp(hemisphere,'CORTEX_LEFT')
    side = 'L';
end

% unix(['~/Applications2/workbench/bin_linux64/wb_command -cifti-separate ' indirRL ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_RL.cortex.' side '.func.gii']);
% unix(['~/Applications2/workbench/bin_linux64/wb_command -cifti-separate ' indirLR ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_LR.cortex.' side '.func.gii']);

% NOTE, if wb_command produces an error, try editing the lines below 
% to contain the full path to wb_command. An example is provided in lines 22 and
% 23
unix(['wb_command -cifti-separate ' indirRL ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_RL.cortex.' side '.func.gii']);
unix(['wb_command -cifti-separate ' indirLR ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_LR.cortex.' side '.func.gii']);

% Example where the full path to wb_command is specified:
%unix(['/home/samsi/Applications2/workbench/bin_linux64/wb_command -cifti-separate ' indirRL ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_RL.cortex.' side '.func.gii']);
%unix(['/home/samsi/Applications2/workbench/bin_linux64/wb_command -cifti-separate ' indirLR ' COLUMN -metric ' hemisphere ' ' tempdir  'temp' '_LR.cortex.' side '.func.gii']);

giiRL = gifti([tempdir 'temp' '_RL.cortex.' side '.func.gii']);

dataRL = giiRL.cdata;

giiLR = gifti([tempdir 'temp' '_LR.cortex.' side '.func.gii']);
dataLR = giiLR.cdata;

unix(['rm ' tempdir  'temp' '_RL.cortex.' side '.func.gii']);
unix(['rm ' tempdir  'temp' '_LR.cortex.' side '.func.gii']);

% Concatenate sessions 1 and 2:
data = [dataRL,dataLR];
load(['./supportingdatafiles/mapGIFTItoCIFTI_cortex_',side,'.mat']);
if strcmp(hemisphere,'CORTEX_RIGHT')
    data = double(data(mapping_CORTEX_RIGHT,:));
else data = double(data(mapping_CORTEX_LEFT,:));
end

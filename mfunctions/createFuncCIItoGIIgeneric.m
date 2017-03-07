function [] = createFuncCIItoGIIgeneric(dmat,hemisphere,filename)
% Benjamin Risk
% This function exports a GII file but does not specify which anatomical
% structure is exported. The structure has to be specified in wb_view.
% Input:
%       dmat 29716 x TR in CORTEX_RIGHT; 29696 in CORTEX_LEFT
%       hemisphere 'CORTEX_RIGHT' or 'CORTEX_LEFT'
%       filename  .func.gii will be appended to "filename"
% Output:
%   filename.func.gii surface file for viewing in wb_view

if nargin==2
    error('createFuncCIItoGIIgeneric requires 3 arguments')
end

[nVertex,nTR] = size(dmat);
if nVertex<nTR
    error('dmat must be nVertex x nTR')
end

if strcmp(hemisphere,'CORTEX_RIGHT')
    load ./supportingdatafiles/mapGIFTItoCIFTI_cortex_R.mat
    shell = zeros(32492,nTR);
    shell(mapping_CORTEX_RIGHT,:) = dmat;
elseif strcmp(hemisphere,'CORTEX_LEFT')
    load ./supportingdatafiles/mapGIFTItoCIFTI_cortex_L.mat
    shell = zeros(32492,nTR);
    shell(mapping_CORTEX_LEFT,:) = dmat;
else 
    error('Implemented for either CORTEX_RIGHT or CORTEX_LEFT only')
end
temp = gifti(shell);
save(temp,[filename '.func.gii'],'Base64Binary');

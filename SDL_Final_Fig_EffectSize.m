
clear
close all

cd ..
SDL.path = pwd; % path of the project
cd('Scripts');
SDL.raw = fullfile(SDL.path,'Original','PGC_Data_Covariates_check_for_delin_withCAPS.xlsx'); % raw data file
SDL.out = fullfile(SDL.path,'Outputs_ComBat_Reduction');
SDL.Boxpath = fullfile('C:\Users\sunde\Documents\MATLAB','SDLToolBox');

fn = fullfile(SDL.out,'CT','Results_EffectSize_CT_corr_.mat');
load(fn,'kk','d'); % kk: the list of the top-n areas, d: effect size of cortical atrophy
load(fullfile(SDL.path,'Original','RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations

k = 148; % top-3 areas
area = [1:k]';
tbl = table(area);
a = cell2mat(D(kk(1:k),2:end));
tbl.x = a(:,1); tbl.y = a(:,2); tbl.z = a(:,3);
tbl.color = sign(d(kk(1:k)))'+2; % PTSD<CONT, 0; PTSD>CONT, 2
tbl.size  = abs(d(kk(1:k)))';
tbl.name  = repmat('-',k,1);
fn = fullfile(SDL.out,'CT','network');
writetable(tbl(:,2:end),[fn,'.dat'],'WriteVariableNames',false,'Delimiter',' ');
movefile([fn,'.dat'],[fn,'.node']); % make the node fle

fn = fullfile(SDL.out,'CT','Data_SC_CT_corr_.mat');
load(fn,'S0');
M = zeros(k,k); % no connections, to show nodes only
% M = S0.Y1(kk(1:k),kk(1:k));
% M(1:1+size(M,1):end) = 0; % convert diagonal into 0s
fn = fullfile(SDL.out,'CT','network');
fid = fopen([fn,'.dat'],'w+');
for ii = 1:size(M,1)
    fprintf(fid,'%1.3f\t',M(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)
movefile([fn,'.dat'],[fn,'.edge']); % make the edge fle
% 
% addpath(fullfile(SDL.Boxpath,'BrainNetViewer_20191031'),'BrainNet');
BrainNet_MapCfg(fullfile(SDL.Boxpath,'BrainNetViewer_20191031','Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv'),...
    [fn,'.node'],[fn,'.edge'],fullfile(SDL.path,'Scripts','FigureOption.mat'));

%BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','network.node','network.edge','FigureOption.mat');

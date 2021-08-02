function SDL_Final_Fig1(SDL)

% make figures for manuscript


%% Fig1A

% =========================Plot CT-based effect size
% 
% --- Nodal size, represents the absoluate value of the effect size of (PTSD versus Non-PTSD) 
% --- Nodal color, red = PTSD > Non-PTSD, blue = PTSD < Non-PTSD
SDL.out1 = fullfile(SDL.path,'Outputs_ComBat_Reduction');
SDL.Boxpath = fullfile('C:\Users\sunde\Documents\MATLAB','SDLToolBox');

fn = fullfile(SDL.out1,'CT','Results_EffectSize_CT_corr_.mat');
load(fn,'kk','d'); % kk: the list of the top-n areas, d: effect size of cortical atrophy
load(fullfile(SDL.path,'Original','RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations

k = 1:148; % top-n areas
area = k';
tbl = table(area);
a = cell2mat(D(:,2:end));
tbl.x = a(:,1); tbl.y = a(:,2); tbl.z = a(:,3);
tbl.color = sign(d(k))'+2; % PTSD<CONT, 0; PTSD>CONT, 2
% tbl.color = ones(length(k),1); % other nodes
% tbl.color(kk(1:10)) = 0; % top-n, PTSD < non-PTSD
% tbl.color(kk(65:84)) = 2;% top-n, PTSD ~ non-PTSD
% tbl.color(kk(139:148)) = 3;% top-n,PTSD > non-PTSD
tbl.size  = abs(d(k))';
tbl.name  = repmat('-',length(k),1);
fn = fullfile(SDL.out,'CT','network');
writetable(tbl(:,2:end),[fn,'.dat'],'WriteVariableNames',false,'Delimiter',' ');
movefile([fn,'.dat'],[fn,'.node']); % make the node fle

fn = fullfile(SDL.out,'CT','Data_SC_CT_corr_.mat');
load(fn,'S0');
M = zeros(length(k),length(k)); % no connections, to show nodes only

% k = 1:3;     % show connectiones of top-3 regions of PTSD < Non-PTSD
% M(k,k) = S0.Y1(kk(k),kk(k));
% k = 146:148; % show connectiones of top-3 regions of PTSD > Non-PTSD
% M(k,k) = S0.Y1(kk(k),kk(k));
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
BrainNet_MapCfg(fullfile(SDL.Boxpath,'BrainNetViewer_20191031','Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv'),...
    [fn,'.node'],[fn,'.edge'],fullfile(SDL.path,'Scripts','FigureOption.mat'));

% ==========================

% ========================= Plot SA-based effect size
% 
% --- Nodal size, represents the absoluate value of the effect size of (PTSD versus Non-PTSD) 
% --- Nodal color, red = PTSD > Non-PTSD, blue = PTSD < Non-PTSD
SDL.out1 = fullfile(SDL.path,'Outputs_ComBat_Reduction');
SDL.Boxpath = fullfile('C:\Users\sunde\Documents\MATLAB','SDLToolBox');

fn = fullfile(SDL.out1,'SA','Results_EffectSize_SA_corr_.mat');
load(fn,'kk','d'); % kk: the list of the top-n areas, d: effect size of cortical atrophy
load(fullfile(SDL.path,'Original','RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations

k = 1:148; % top-3 areas
area = k';
tbl = table(area);
a = cell2mat(D(:,2:end));
tbl.x = a(:,1); tbl.y = a(:,2); tbl.z = a(:,3);
tbl.color = sign(d(k))'+2; % PTSD<CONT, 0; PTSD>CONT, 2
% tbl.color = ones(length(k),1); % other nodes
% tbl.color(kk(1:10)) = 0; % top-n, PTSD < non-PTSD
% tbl.color(kk(65:84)) = 2;% top-n, PTSD ~ non-PTSD
% tbl.color(kk(139:148)) = 3;% top-n,PTSD > non-PTSD
tbl.size  = abs(d(k))';
tbl.name  = repmat('-',length(k),1);
fn = fullfile(SDL.out,'SA','network');
writetable(tbl(:,2:end),[fn,'.dat'],'WriteVariableNames',false,'Delimiter',' ');
movefile([fn,'.dat'],[fn,'.node']); % make the node fle

fn = fullfile(SDL.out,'SA','Data_SC_SA_corr_.mat');
load(fn,'S0');
M = zeros(length(k),length(k)); % no connections, to show nodes only

% k = 1:3;     % show connectiones of top-3 regions of PTSD < Non-PTSD
% M(k,k) = S0.Y1(kk(k),kk(k));
% k = 146:148; % show connectiones of top-3 regions of PTSD > Non-PTSD
% M(k,k) = S0.Y1(kk(k),kk(k));
% M(1:1+size(M,1):end) = 0; % convert diagonal into 0s
fn = fullfile(SDL.out,'SA','network');
fid = fopen([fn,'.dat'],'w+');
for ii = 1:size(M,1)
    fprintf(fid,'%1.3f\t',M(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)
movefile([fn,'.dat'],[fn,'.edge']); % make the edge fle
% 
BrainNet_MapCfg(fullfile(SDL.Boxpath,'BrainNetViewer_20191031','Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv'),...
    [fn,'.node'],[fn,'.edge'],fullfile(SDL.path,'Scripts','FigureOption.mat'));

% % ==========================





% Upper panel
Ana = {'CT', 'corr',        'und', {},{},{},''};
SDL.data_type = Ana(1,1); % CT or SA
SDL.ana_type  = Ana(1,2); % corr, partialcorr or med
SDL.XYM       = Ana(1,4:7); % for mediation analyses only

fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); % index of two groups
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'S0');
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    X1(i-1) = SC0(i).mean1;
    X2(i-1) = SC0(i).mean2;
    p1(i-1) = SC1(i).p1;
    p2(i-1) = SC1(i).p2;
    
    X1R(i-1) = mean(SC1(i).mean1);
    X2R(i-1) = mean(SC1(i).mean2);
    
    X1RA(i-1,:) = SC1(i).mean1;
    X2RA(i-1,:) = SC1(i).mean2;

    CI1L(i-1) = SC1(i).CI1(1);
    CI1R(i-1) = SC1(i).CI1(2);
    
    CI2L(i-1) = SC1(i).CI2(1);
    CI2R(i-1) = SC1(i).CI2(2);
    
    diffL(i-1) = SC1(i).diffCI(1);
    diffR(i-1) = SC1(i).diffCI(2);
    p_diff(i-1) = SC1(i).p_diff;
end

% get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure; % original vs permutation
subplot(2,2,1);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [CI1L, fliplr(CI1R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,X1,'r-',X,X1R,'b-','LineWidth',2);
% xlabel('Top-N Regions'); ylabel('CT-Based Structural Covariance');%legend('95% CI','Areas with Reduction','Random Areas'); 
% title('PTSD');
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',12); set(gca,'XTick',0:25:150);
ylabel('Mean SC','FontName', 'Arial','FontSize',12);%CT-bsed
title('PTSD','FontName', 'Arial','FontSize',12);
box off

subplot(2,2,2);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [CI2L, fliplr(CI2R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,X2,'r-',X,X2R,'b-','LineWidth',2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',12); set(gca,'XTick',0:25:150);
ylabel('Mean SC','FontName', 'Arial','FontSize',12);% CT-based
title('CONT','FontName', 'Arial','FontSize',12);
box off


% Lower panel
Ana = {'SA', 'corr',        'und', {},{},{},''};
SDL.data_type = Ana(1,1); % CT or SA
SDL.ana_type  = Ana(1,2); % corr, partialcorr or med
SDL.XYM       = Ana(1,4:7); % for mediation analyses only

fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); % index of two groups
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'S0');
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    X1(i-1) = SC0(i).mean1;
    X2(i-1) = SC0(i).mean2;
    p1(i-1) = SC1(i).p1;
    p2(i-1) = SC1(i).p2;
    
    X1R(i-1) = mean(SC1(i).mean1);
    X2R(i-1) = mean(SC1(i).mean2);
    
    X1RA(i-1,:) = SC1(i).mean1;
    X2RA(i-1,:) = SC1(i).mean2;

    CI1L(i-1) = SC1(i).CI1(1);
    CI1R(i-1) = SC1(i).CI1(2);
    
    CI2L(i-1) = SC1(i).CI2(1);
    CI2R(i-1) = SC1(i).CI2(2);
    
    diffL(i-1) = SC1(i).diffCI(1);
    diffR(i-1) = SC1(i).diffCI(2);
    p_diff(i-1) = SC1(i).p_diff;
end

subplot(2,2,3);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [CI1L, fliplr(CI1R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,X1,'r-',X,X1R,'b-','LineWidth',2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',12); set(gca,'XTick',0:25:150);
ylabel('Mean SC','FontName', 'Arial','FontSize',12); % SA-based network
legend('95% CI','Real Network','Random Network'); 
title('PTSD','FontName', 'Arial','FontSize',12);
box off

subplot(2,2,4);
X = 2:size(S0.Y1,1); x2 = [X, fliplr(X)]; inBetween = [CI2L, fliplr(CI2R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,X2,'r-',X,X2R,'b-','LineWidth',2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',12); set(gca,'XTick',0:25:150);
ylabel('Mean SC','FontName', 'Arial','FontSize',12); % SA-based network
title('CONT','FontName', 'Arial','FontSize',12);
box off

savefig(fullfile(SDL.out,'Figure 1.fig'));



% brain map - left
k = 50; % top-3 areas
fn = fullfile(SDL.out,'CT','Results_EffectSize_CT_corr_.mat');
load(fn,'kk','d'); % kk: the list of the top-n areas, d: effect size of cortical atrophy
load(fullfile(SDL.path,'Original','RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations

area = [1:148]';
tbl = table(area);
a = cell2mat(D(:,2:end));
tbl.x = a(:,1); tbl.y = a(:,2); tbl.z = a(:,3);
tbl.color = ones(148,1);
tbl.size  = d';
tbl.name  = repmat('-',148,1);
tbl = tbl(kk(1:k),:);
fn = fullfile(SDL.out,'CT','network');
writetable(tbl(:,2:end),[fn,'.dat'],'WriteVariableNames',false,'Delimiter',' ');
movefile([fn,'.dat'],[fn,'.node']); % make the node fle

fn = fullfile(SDL.out,'CT','Data_SC_CT_corr_.mat');
load(fn,'S0');
M = S0.Y1(kk(1:k),kk(1:k));
M(1:1+size(M,1):end) = 0; % convert diagonal into 0s
fn = fullfile(SDL.out,'CT','network');
fid = fopen([fn,'.dat'],'w+');
for ii = 1:size(M,1)
    fprintf(fid,'%1.3f\t',M(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)
movefile([fn,'.dat'],[fn,'.edge']); % make the edge fle

addpath(fullfile(SDL.Boxpath,'BrainNetViewer_20191031'),'BrainNet');
BrainNet_MapCfg(fullfile(SDL.Boxpath,'BrainNetViewer_20191031','Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv'),...
    [fn,'.node'],[fn,'.edge'],fullfile(SDL.path,'Scripts','FigureOption.mat'));




% brain map - right
fn = fullfile(SDL.out,'SA','Results_EffectSize_SA_corr_.mat');
load(fn,'kk','d'); % kk: the list of the top-n areas, d: effect size of cortical atrophy
load(fullfile(SDL.path,'Original','RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations

k = 13; % top-3 areas for SA decrease in PTSD vs CONT; but top-13 areas for SA increase
area = [1:k]';
tbl = table(area);
a = cell2mat(D(kk(1:k),2:end));
tbl.x = a(:,1); tbl.y = a(:,2); tbl.z = a(:,3);
tbl.color = ones(k,1);
tbl.size(kk(k))  = abs(d(kk(k)))';
tbl.name  = repmat('name',k,1);
fn = fullfile(SDL.out,'SA','network');
writetable(tbl(:,2:end),[fn,'.dat'],'WriteVariableNames',false,'Delimiter',' ');
movefile([fn,'.dat'],[fn,'.node']); % make the node fle

fn = fullfile(SDL.out,'SA','Data_SC_SA_corr_.mat');
load(fn,'S0');
M = S0.Y1(kk(1:k),kk(1:k));
M(1:1+size(M,1):end) = 0; % convert diagonal into 0s
fn = fullfile(SDL.out,'SA','network');
fid = fopen([fn,'.dat'],'w+');
for ii = 1:size(M,1)
    fprintf(fid,'%1.3f\t',M(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)
movefile([fn,'.dat'],[fn,'.edge']); % make the edge fle

addpath(fullfile(SDL.Boxpath,'BrainNetViewer_20191031'),'BrainNet');
BrainNet_MapCfg(fullfile(SDL.Boxpath,'BrainNetViewer_20191031','Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv'),...
    [fn,'.node'],[fn,'.edge'],fullfile(SDL.path,'Scripts','FigureOption.mat'));




%% End
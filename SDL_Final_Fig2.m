function SDL_Final_Fig2(SDL)

% make figures for manuscript



% Left panel
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


figure;
a = [1:2]';
b = sum([(X1-X1R)',...
    (X2-X2R)'])';
b= bar(0,b);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',16)
% xlabel('PTSD Severity','FontName', 'Arial','FontSize',14); 
ylabel('AUC','FontName', 'Arial','FontSize',16);
xticklabels({' '});
box off
b(1).FaceColor = 'r';
b(2).FaceColor = 'g';




figure; % original vs permutation
% subplot(1,2,1);
X = 2:size(S0.Y1,1); 
% x2 = [X, fliplr(X)]; inBetween = [diffL, fliplr(diffR)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,(X1-X1R),'r-',X,(X2-X2R),'g-','LineWidth',2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',14)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',18); set(gca,'XTick',0:10:150);
ylabel('CT-Based Structural Covariance','FontName', 'Arial','FontSize',18);%legend('95% CI','Areas with Reduction','Random Areas'); 
title('Real-Random','FontName', 'Arial','FontSize',18);
legend('PTSD','CONT','FontName', 'Arial','FontSize',16);
box off




% brain map - left
fn = fullfile(SDL.out,'CT','Results_EffectSize_CT_corr_.mat');
load(fn,'kk','d'); % kk: the list of the top-n areas, d: effect size of cortical atrophy
load(fullfile(SDL.path,'Original','RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations

k = 46; % top-31 or 78 areas
area = [1:k]';
tbl = table(area);
a = cell2mat(D(kk(1:k),2:end));
tbl.x = a(:,1); tbl.y = a(:,2); tbl.z = a(:,3);
tbl.color = ones(k,1);
tbl.size  = abs(d(kk(1:k)))';
tbl.name  = repmat('name',k,1);
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

BoxPath = '/Users/ds366/Documents/MATLAB/SDLToolBox';
addpath(fullfile(BoxPath,'BrainNetViewer_20171031'),'BrainNet');
BrainNet_MapCfg('/Users/ds366/Documents/MATLAB/SDLToolBox/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
    [fn,'.node'],[fn,'.edge'],fullfile(SDL.path,'Scripts','FigureOption.mat'));













% Right panel
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

figure;
a = [1:2]';
b = sum([(X1-X1R)',...
    (X2-X2R)'])';
b= bar(0,b);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',16)
% xlabel('PTSD Severity','FontName', 'Arial','FontSize',14); 
ylabel('AUC','FontName', 'Arial','FontSize',16);
xticklabels({' '});
box off
b(1).FaceColor = 'r';
b(2).FaceColor = 'g';




figure; % original vs permutation
% subplot(1,2,1);
X = 2:size(S0.Y1,1); 
% x2 = [X, fliplr(X)]; inBetween = [diffL, fliplr(diffR)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
plot(X,(X1-X1R),'r-',X,(X2-X2R),'g-','LineWidth',2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',14)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',18); set(gca,'XTick',0:10:150);
ylabel('SA-Based Structural Covariance','FontName', 'Arial','FontSize',18);%legend('95% CI','Areas with Reduction','Random Areas'); 
title('Real-Random','FontName', 'Arial','FontSize',18);
legend('PTSD','CONT','FontName', 'Arial','FontSize',16);
box off

savefig(fullfile(SDL.out,'Figure 2.fig'));


%% End
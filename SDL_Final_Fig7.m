function SDL_Final_Fig7(SDL)

% make figures for manuscript



% Left panel
Ana = {'CT', 'corr',        'und', {},{},{},''};
i=1;
SDL.data_type = Ana(i,1); % CT or SA
SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
SDL.XYM       = Ana(i,4:7); % for mediation analyses only
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); % index of two groups

% Effect size
fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn,'kk','d'); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);

S0 = SDL_SCA(R,idx1,idx2,@corr,0);     % SC of original PTSD & CONT data


% Only pos SC
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
% figure('DefaultAxesFontSize',18)
figure();
ax1 = subplot(2,2,1); % PTSD
a = triu(S0.Y1,1); a(a<0) = NaN; % only pos SC
a = a + a'; % structural covariance after removing Inf in the diagnose
a = [d',nanmean(a)']; 
%a = a(a(:,1)<0,:);
scatter(ax1,a(:,1),a(:,2),45,'filled');h1 = lsline(ax1); h1.Color = 'r'; h1.LineWidth = 2; %axis([-0.11,0.01,-0.025,0.025]);
title('PTSD'); xlabel('Effect Size'); ylabel('CT-Based Mean SC');
[r,p] = corr(a(:,1),a(:,2));
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(-0.1,0.022,txt,'FontSize',18);

ax2 = subplot(2,2,3); % CONT
a = triu(S0.Y2,1); a(a<0) = NaN; % only pos SC
a = a + a'; % structural covariance after removing Inf in the diagnose
a = [d',nanmean(a)']; 
%a = a(a(:,1)<0,:);
scatter(ax2,a(:,1),a(:,2),45,'filled');h2 = lsline(ax2); h2.Color = 'r'; h2.LineWidth = 2;  %axis([-0.11,0.01,-0.025,0.025]);
title('CONT'); xlabel('Effect Size'); ylabel('CT-Based Structural Covariance');
[r,p] = corr(a(:,1),a(:,2));
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(-0.1,0.022,txt,'FontSize',18);


% Right panel
Ana = {'SA', 'corr',        'und', {},{},{},''};
i=1;
SDL.data_type = Ana(i,1); % CT or SA
SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
SDL.XYM       = Ana(i,4:7); % for mediation analyses only
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); % index of two groups

% Effect size
fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn,'kk','d'); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);

S0 = SDL_SCA(R,idx1,idx2,@corr,0);     % SC of original PTSD & CONT data


% Only neg efect size

ax1 = subplot(2,2,2); % PTSD
a = triu(S0.Y1,1); a(a<0) = NaN; % only pos SC
a = a + a'; % structural covariance after removing Inf in the diagnose
a = [d',nanmean(a)']; 
%a = a(a(:,1)<0,:);
scatter(ax1,a(:,1),a(:,2),45,'filled');h1 = lsline(ax1); h1.Color = 'r'; h1.LineWidth = 2; %axis([-0.11,0.01,-0.025,0.025]);
title('PTSD'); xlabel('Effect Size'); ylabel('SA-Based Mean SC');
[r,p] = corr(a(:,1),a(:,2));
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(-0.1,0.022,txt,'FontSize',18);

ax2 = subplot(2,2,4); % CONT
a = triu(S0.Y2,1); a(a<0) = NaN; % only pos SC
a = a + a'; % structural covariance after removing Inf in the diagnose
a = [d',nanmean(a)']; 
%a = a(a(:,1)<0,:);
scatter(ax2,a(:,1),a(:,2),45,'filled');h2 = lsline(ax2); h2.Color = 'r'; h2.LineWidth = 2;  %axis([-0.11,0.01,-0.025,0.025]);
title('CONT'); xlabel('Effect Size'); ylabel('SA-Based Mean');
[r,p] = corr(a(:,1),a(:,2));
txt = sprintf('R=%1.3f, p=%1.3f',r,p);
text(-0.1,0.022,txt,'FontSize',18);


%% End
% savefig(fullfile(SDL.out,'Figure 7.fig'));
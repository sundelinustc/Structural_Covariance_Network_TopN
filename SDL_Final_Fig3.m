function SDL_Final_Fig3(SDL)

% make figures for manuscript


%% Figure 3. PTSD severity

% Left panel
Ana = {'CT_PTSDsev10', 'corr',    'und', {},{},{},''}; i=1;
SDL.data_type = Ana(i,1); % CT or SA
SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
SDL.XYM       = Ana(i,4:7); % for mediation analyses only
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Results_TopN_SC_CI_p_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    
    X1(i-1) = SC0(i).mean1; 
    X2(i-1) = SC0(i).mean2; 
    X3(i-1) = SC0(i).mean3; 
    X4(i-1) = SC0(i).mean4; 
    X5(i-1) = SC0(i).mean5; 
    
    p1(i-1) = SC1(i).p1; 
    p2(i-1) = SC1(i).p2; 
    p3(i-1) = SC1(i).p3; 
    p4(i-1) = SC1(i).p4;
    p5(i-1) = SC1(i).p5; 
    
    X1R(i-1) = mean(SC1(i).mean1); 
    X2R(i-1) = mean(SC1(i).mean2); 
    X3R(i-1) = mean(SC1(i).mean3);
    X4R(i-1) = mean(SC1(i).mean4);
    X5R(i-1) = mean(SC1(i).mean5);
    
    
    
    X1RA(i-1,:) = SC1(i).mean1; 
    X2RA(i-1,:) = SC1(i).mean2; 
    X3RA(i-1,:) = SC1(i).mean3; 
    X4RA(i-1,:) = SC1(i).mean4; 
    X5RA(i-1,:) = SC1(i).mean5; 
end

% get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure; % original vs permutation

subplot(2,2,1);
L = 1:length(X1);
plot(L,(X1-X1R),...
    L,(X2-X2R),...
    L,(X3-X3R),...
    L,(X4-X4R),...
    L,(X5-X5R),'LineWidth', 2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
% legend('  0<=CAPS<2','  2<=CAPS<18','18<=CAPS<44','44<=CAPS<67','67<=CAPS<120');
xlabel('Top-N Regions','FontName', 'Arial','FontSize',14); set(gca,'XTick',0:25:150);
ylabel('CT-Based Structural Covariance','FontName', 'Arial','FontSize',14);
title('Real - Random','FontName', 'Arial','FontSize',14)
box off


subplot(2,2,3);
a = [1:5]';
b = mean([(X1-X1R)',...
    (X2-X2R)',...
    (X3-X3R)',...
    (X4-X4R)',...
    (X5-X5R)'])';
bar(0,b);
% hold on
% f = fit(a,b,'poly23');
% plot(f,a,b);
xlabel('PTSD Severity','FontName', 'Arial','FontSize',14); 
ylabel('AUC','FontName', 'Arial','FontSize',14);
xticklabels({' '});
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',1)
% set(gca,'XTickLabelMode','auto')
box off












% Right panel
Ana = {'SA_PTSDsev10', 'corr',    'und', {},{},{},''}; i=1;
SDL.data_type = Ana(i,1); % CT or SA
SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
SDL.XYM       = Ana(i,4:7); % for mediation analyses only
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Results_TopN_SC_CI_p_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    
    X1(i-1) = SC0(i).mean1; 
    X2(i-1) = SC0(i).mean2; 
    X3(i-1) = SC0(i).mean3; 
    X4(i-1) = SC0(i).mean4; 
    X5(i-1) = SC0(i).mean5; 
    
    p1(i-1) = SC1(i).p1; 
    p2(i-1) = SC1(i).p2; 
    p3(i-1) = SC1(i).p3; 
    p4(i-1) = SC1(i).p4;
    p5(i-1) = SC1(i).p5; 
    
    X1R(i-1) = mean(SC1(i).mean1); 
    X2R(i-1) = mean(SC1(i).mean2); 
    X3R(i-1) = mean(SC1(i).mean3);
    X4R(i-1) = mean(SC1(i).mean4);
    X5R(i-1) = mean(SC1(i).mean5);
    
    
    
    X1RA(i-1,:) = SC1(i).mean1; 
    X2RA(i-1,:) = SC1(i).mean2; 
    X3RA(i-1,:) = SC1(i).mean3; 
    X4RA(i-1,:) = SC1(i).mean4; 
    X5RA(i-1,:) = SC1(i).mean5; 
end

% get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
% figure; % original vs permutation

subplot(2,2,2);
L = 1:length(X1);
plot(L,(X1-X1R),...
    L,(X2-X2R),...
    L,(X3-X3R),...
    L,(X4-X4R),...
    L,(X5-X5R),'LineWidth', 2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
legend('  0<=CAPS<4','  4<=CAPS<20','20<=CAPS<46','46<=CAPS<68','68<=CAPS<120');
% legend('  0<=CAPS<2','  2<=CAPS<18','18<=CAPS<44','44<=CAPS<67','67<=CAPS<120');
xlabel('Top-N Regions','FontName', 'Arial','FontSize',14); set(gca,'XTick',0:25:150);
ylabel('SA-Based Structural Covariance','FontName', 'Arial','FontSize',14);
title('Real - Random','FontName', 'Arial','FontSize',14)
box off


subplot(2,2,4);
a = [1:5]';
b = mean([(X1-X1R)',...
    (X2-X2R)',...
    (X3-X3R)',...
    (X4-X4R)',...
    (X5-X5R)'])';
bar(0,b);
% hold on
% f = fit(a,b,'poly23');
% plot(f,a,b);
xlabel('PTSD Severity','FontName', 'Arial','FontSize',14); 
ylabel('AUC','FontName', 'Arial','FontSize',14);
xticklabels({' '});
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',1)
% set(gca,'XTickLabelMode','auto')
box off



%% End
savefig(fullfile(SDL.out,'Figure 3.fig'));
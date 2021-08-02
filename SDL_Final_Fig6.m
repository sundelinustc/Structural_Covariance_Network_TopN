function SDL_Final_Fig6(SDL)

% make figures for manuscript


%% Figure 6. PTSD diagnosis x Age

% Left panel
Ana = {'CT_Age10', 'corr',    'und', {},{},{},''}; i=1;
SDL.data_type = Ana(i,1); % CT or SA
SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
SDL.XYM       = Ana(i,4:7); % for mediation analyses only
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    
X11(i-1) = SC0(i).mean11; 
    X12(i-1) = SC0(i).mean12; 
    X13(i-1) = SC0(i).mean13; 
    X14(i-1) = SC0(i).mean14; 
    X15(i-1) = SC0(i).mean15; 
    X16(i-1) = SC0(i).mean16; 
    X17(i-1) = SC0(i).mean17;
    X18(i-1) = SC0(i).mean18;
    
    X21(i-1) = SC0(i).mean21; 
    X22(i-1) = SC0(i).mean22;
    X23(i-1) = SC0(i).mean23;
    X24(i-1) = SC0(i).mean24;
    X25(i-1) = SC0(i).mean25;
    X26(i-1) = SC0(i).mean26;
    X27(i-1) = SC0(i).mean27;
    X28(i-1) = SC0(i).mean28;
        
    X11R(i-1) = mean(SC1(i).mean11); 
    X12R(i-1) = mean(SC1(i).mean12); 
    X13R(i-1) = mean(SC1(i).mean13);
    X14R(i-1) = mean(SC1(i).mean14);
    X15R(i-1) = mean(SC1(i).mean15);
    X16R(i-1) = mean(SC1(i).mean16);
    X17R(i-1) = mean(SC1(i).mean17);
    X18R(i-1) = mean(SC1(i).mean18);
    
    X21R(i-1) = mean(SC1(i).mean21); 
    X22R(i-1) = mean(SC1(i).mean22);
    X23R(i-1) = mean(SC1(i).mean23);
    X24R(i-1) = mean(SC1(i).mean24);
    X25R(i-1) = mean(SC1(i).mean25);
    X26R(i-1) = mean(SC1(i).mean26);
    X27R(i-1) = mean(SC1(i).mean27);
    X28R(i-1) = mean(SC1(i).mean28);
end


% subplot(2,2,3);
figure;
a = [1:8]';
b1 = sum([(X11-X11R)',...
    (X12-X12R)',...
    (X13-X13R)',...
    (X14-X14R)',...
    (X15-X15R)',...
    (X16-X16R)',...
    (X17-X17R)',...
    (X18-X18R)'])';
b2 = sum([(X21-X21R)',...
    (X22-X22R)',...
    (X23-X23R)',...
    (X24-X24R)',...
    (X25-X25R)',...
    (X26-X26R)',...
    (X27-X27R)',...
    (X28-X28R)'])';
b= bar([b1,b2]);
b(1).FaceColor = 'r'; % PTSD
b(2).FaceColor = 'g'; % CONT
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',14)
xlabel('Age','FontName', 'Arial','FontSize',16); 
ylabel('AUC','FontName', 'Arial','FontSize',16);
xticklabels({'<10','10~15','15~20','20~30','30~40','40~50','50~60','>=60'});
legend('PTSD','Controls','FontName', 'Arial','FontSize',16);
box off


figure; % original vs permutation
% subplot(2,2,1);
L = 1:length(X11);
plot(L,(X11-X11R)-(X21-X21R),...
    L,(X12-X12R)-(X22-X22R),...
    L,(X13-X13R)-(X23-X23R),...
    L,(X14-X14R)-(X24-X24R),...
    L,(X15-X15R)-(X26-X26R),...
    L,(X16-X16R)-(X26-X26R),...
    L,(X17-X17R)-(X27-X27R),...
    L,(X18-X18R)-(X28-X28R),'LineWidth', 2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
ylim([-0.15,0.16]);
% legend('<20','>=20 and <30','>=30 and <40','>=40 and <50','>=50 and <60','>=60');
xlabel('Top-N Regions','FontName', 'Arial','FontSize',12);  set(gca,'XTick',0:10:150);
ylabel('CT-Based Structural Covariance','FontName', 'Arial','FontSize',12);
title('(PTSD-CONT)_R_e_a_l_-_R_a_n_d_o_m','FontName', 'Arial','FontSize',12);
legend('  0<=Age<10','10<=Age<15','15<=Age<20','20<=Age<30','30<=Age<40','40<=Age<50','50<=Age<60','60<=Age');
box off





% Right panel
Ana = {'SA_Age10', 'corr',    'und', {},{},{},''}; i=1;
SDL.data_type = Ana(i,1); % CT or SA
SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
SDL.XYM       = Ana(i,4:7); % for mediation analyses only
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    
X11(i-1) = SC0(i).mean11; 
    X12(i-1) = SC0(i).mean12; 
    X13(i-1) = SC0(i).mean13; 
    X14(i-1) = SC0(i).mean14; 
    X15(i-1) = SC0(i).mean15; 
    X16(i-1) = SC0(i).mean16; 
    X17(i-1) = SC0(i).mean17;
    X18(i-1) = SC0(i).mean18;
    
    X21(i-1) = SC0(i).mean21; 
    X22(i-1) = SC0(i).mean22;
    X23(i-1) = SC0(i).mean23;
    X24(i-1) = SC0(i).mean24;
    X25(i-1) = SC0(i).mean25;
    X26(i-1) = SC0(i).mean26;
    X27(i-1) = SC0(i).mean27;
    X28(i-1) = SC0(i).mean28;
        
    X11R(i-1) = mean(SC1(i).mean11); 
    X12R(i-1) = mean(SC1(i).mean12); 
    X13R(i-1) = mean(SC1(i).mean13);
    X14R(i-1) = mean(SC1(i).mean14);
    X15R(i-1) = mean(SC1(i).mean15);
    X16R(i-1) = mean(SC1(i).mean16);
    X17R(i-1) = mean(SC1(i).mean17);
    X18R(i-1) = mean(SC1(i).mean18);
    
    X21R(i-1) = mean(SC1(i).mean21); 
    X22R(i-1) = mean(SC1(i).mean22);
    X23R(i-1) = mean(SC1(i).mean23);
    X24R(i-1) = mean(SC1(i).mean24);
    X25R(i-1) = mean(SC1(i).mean25);
    X26R(i-1) = mean(SC1(i).mean26);
    X27R(i-1) = mean(SC1(i).mean27);
    X28R(i-1) = mean(SC1(i).mean28);
end

% subplot(2,2,3);
figure;
a = [1:6]';
b1 = sum([(X11-X11R)',...
    (X12-X12R)',...
    (X13-X13R)',...
    (X14-X14R)',...
    (X15-X15R)',...
    (X16-X16R)',...
    (X17-X17R)',...
    (X18-X18R)'])';
b2 = sum([(X21-X21R)',...
    (X22-X22R)',...
    (X23-X23R)',...
    (X24-X24R)',...
    (X25-X25R)',...
    (X26-X26R)',...
    (X27-X27R)',...
    (X28-X28R)'])';
b= bar([b1,b2]);
b(1).FaceColor = 'r'; % PTSD
b(2).FaceColor = 'g'; % CONT
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',14)
xlabel('Age','FontName', 'Arial','FontSize',16); 
ylabel('AUC','FontName', 'Arial','FontSize',16);
xticklabels({'<10','10~15','15~20','20~30','30~40','40~50','50~60','>=60'});
legend('PTSD','Controls','FontName', 'Arial','FontSize',16);
box off


figure; % original vs permutation
% subplot(2,2,1);
L = 1:length(X11);
plot(L,(X11-X11R)-(X21-X21R),...
    L,(X12-X12R)-(X22-X22R),...
    L,(X13-X13R)-(X23-X23R),...
    L,(X14-X14R)-(X24-X24R),...
    L,(X15-X15R)-(X26-X26R),...
    L,(X16-X16R)-(X26-X26R),...
    L,(X17-X17R)-(X27-X27R),...
    L,(X18-X18R)-(X28-X28R),'LineWidth', 2);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',12)
ylim([-0.3,0.6]);
% legend('<20','>=20 and <30','>=30 and <40','>=40 and <50','>=50 and <60','>=60');
xlabel('Top-N Regions','FontName', 'Arial','FontSize',12);  set(gca,'XTick',0:10:150);
ylabel('SA-Based Structural Covariance','FontName', 'Arial','FontSize',12);
title('(PTSD-CONT)_R_e_a_l_-_R_a_n_d_o_m','FontName', 'Arial','FontSize',12);
legend('  0<=Age<10','10<=Age<15','15<=Age<20','20<=Age<30','30<=Age<40','40<=Age<50','50<=Age<60','60<=Age');
box off



%% End
savefig(fullfile(SDL.out,'Figure 6.fig'));
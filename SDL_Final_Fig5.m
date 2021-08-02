function SDL_Final_Fig5(SDL)

% make figures for manuscript


%% Figure 3. PTSD diagnosis x Sex

% left panel
Ana = {'CT_Gender', 'corr', 'und', {},{},{},''};
SDL.data_type = Ana(1,1); % CT or SA
SDL.ana_type  = Ana(1,2); % corr, partialcorr or med
SDL.XYM       = Ana(1,4:7); % for mediation analyses only

fdir = fullfile(SDL.out,SDL.data_type{1}(1:2)); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
if strcmp(SDL.data_type{1}(4:end),'Age')
    idx11 = strcmp(T.Group,'PTSD') & (T.Age<nanmedian(T.Age));  % group 11, PTSD & Age < 33 (median)
    idx12 = strcmp(T.Group,'PTSD') & (T.Age>=nanmedian(T.Age)); % group 12, PTSD & Age >= 33 (median)
    idx21 = strcmp(T.Group,'CONT') & (T.Age<nanmedian(T.Age));  % group 21, CONT & Age < 33 (median)
    idx22 = strcmp(T.Group,'CONT') & (T.Age>=nanmedian(T.Age)); % group 22, CONT & Age >= 33 (median)
elseif strcmp(SDL.data_type{1}(4:end),'Gender')
    idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'M');  % group 11, PTSD & Male
    idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'F');  % group 12, PTSD & Female
    idx21 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'M');  % group 21, CONT & Male
    idx22 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'F');  % group 22, CONT & Female
elseif strcmp(SDL.data_type{1}(4:end),'Dep')
%     try
%         idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'1');  % group 11, PTSD & depression
%         idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'0');  % group 12, PTSD & no-depression
%         idx21 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'1');  % group 21, CONT & depression
%         idx22 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'0');  % group 22, CONT & no-depression
%     catch
        idx11 = strcmp(T.Group,'PTSD') & (T.Dep==1);  % group 11, PTSD & depression
        idx12 = strcmp(T.Group,'PTSD') & (T.Dep==0);  % group 12, PTSD & no-depression
        idx21 = strcmp(T.Group,'CONT') & (T.Dep==1);  % group 21, CONT & depression
        idx22 = strcmp(T.Group,'CONT') & (T.Dep==0);  % group 22, CONT & no-depression
%     end
else
end
if strcmp(SDL.ana_type{1},'corr')
    func = @corr;
elseif strcmp(SDL.ana_type{1},'partialcorr')
    func = @partialcorr;
else
end
S0 = SDL_SCA4(R,idx11,idx12,idx21,idx22,func,0);     % SC of original PTSD & CONT data
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    X11(i-1) = SC0(i).mean11; X12(i-1) = SC0(i).mean12; X21(i-1) = SC0(i).mean21; X22(i-1) = SC0(i).mean22;
    
    p11(i-1) = SC1(i).p11; p12(i-1) = SC1(i).p12; p21(i-1) = SC1(i).p21; p22(i-1) = SC1(i).p22;
    
    X11R(i-1) = mean(SC1(i).mean11); X12R(i-1) = mean(SC1(i).mean12); X21R(i-1) = mean(SC1(i).mean21); X22R(i-1) = mean(SC1(i).mean22);
    X11RA(i-1,:) = SC1(i).mean11; X12RA(i-1,:) = SC1(i).mean12; X21RA(i-1,:) = SC1(i).mean21; X22RA(i-1,:) = SC1(i).mean22;
    
    CI11L(i-1) = SC1(i).CI11(1); CI11R(i-1) = SC1(i).CI11(2); CI12L(i-1) = SC1(i).CI12(1); CI12R(i-1) = SC1(i).CI12(2);
    CI21L(i-1) = SC1(i).CI21(1); CI21R(i-1) = SC1(i).CI21(2); CI22L(i-1) = SC1(i).CI22(1); CI22R(i-1) = SC1(i).CI22(2);

    diff_11v21L(i-1) = SC1(i).diffCI_11v21(1); diff_11v21R(i-1) = SC1(i).diffCI_11v21(2);
    diff_12v22L(i-1) = SC1(i).diffCI_12v22(1); diff_12v22R(i-1) = SC1(i).diffCI_12v22(2);
    diff_11v12L(i-1) = SC1(i).diffCI_11v12(1); diff_11v12R(i-1) = SC1(i).diffCI_11v12(2);
    diff_21v22L(i-1) = SC1(i).diffCI_21v22(1); diff_21v22R(i-1) = SC1(i).diffCI_21v22(2);
    diff_main1L(i-1)  = SC1(i).diffCI_main1(1); diff_main1R(i-1) = SC1(i).diffCI_main1(2);
    diff_main2L(i-1)  = SC1(i).diffCI_main2(1); diff_main2R(i-1) = SC1(i).diffCI_main2(2);
    diff_interL(i-1)  = SC1(i).diffCI_inter(1); diff_interR(i-1) = SC1(i).diffCI_inter(2);
    
    
    p_diff_11v21(i-1) = SC1(i).p_diff_11v21; 
    p_diff_12v22(i-1) = SC1(i).p_diff_12v22;
    p_diff_11v12(i-1) = SC1(i).p_diff_11v12;
    p_diff_21v22(i-1) = SC1(i).p_diff_21v22;
    p_diff_main1(i-1) = SC1(i).p_diff_main1;
    p_diff_main2(i-1) = SC1(i).p_diff_main2;
    p_diff_inter(i-1) = SC1(i).p_diff_inter;
end


figure;
a = [1:4]';
b = sum([(X11-X11R)',...
    (X12-X12R)',...
    (X21-X21R)',...
    (X22-X22R)'])';
b= bar(0,b);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',16)
% xlabel('PTSD Severity','FontName', 'Arial','FontSize',14); 
ylabel('AUC','FontName', 'Arial','FontSize',16);
xticklabels({' '});
box off
b(1).FaceColor = 'r';
b(2).FaceColor = 'm';
b(3).FaceColor = 'b';
b(4).FaceColor = 'g';


figure;
% txt1 = 'Dep'; txt2 = 'NoDep';
% subplot(1,2,2);
X = 2:size(S0.Y11,1); 
% x2 = [X, fliplr(X)]; inBetween = [diff_interL, fliplr(diff_interR)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11-X12-X21+X22),'r-',X,(X11R-X12R-X21R+X22R),'b-','LineWidth',2); % interaction
plot(X,(X11-X11R),'r',X,(X12-X12R),'m',X,(X21-X21R),'b',X,(X22-X22R),'g','LineWidth',2); % four groups
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',14)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',16); 
set(gca,'XTick',0:10:150); 
% set(gca,'XTick',0:1:20); 
ylabel('CT-Based Structural Covariance','FontName', 'Arial','FontSize',16);%legend('95% CI','Areas with Reduction','Random Areas'); 
title('Real-Random','FontName', 'Arial','FontSize',16);
% legend('95% CI','Real Network','Random Network','FontName', 'Arial','FontSize',16);
legend('PTSD & Male','PTSD & Female','CONT & Male','CONT & Female','FontName', 'Arial','FontSize',16);
box off








% right panel
Ana = {'SA_Gender', 'corr', 'und', {},{},{},''};
SDL.data_type = Ana(1,1); % CT or SA
SDL.ana_type  = Ana(1,2); % corr, partialcorr or med
SDL.XYM       = Ana(1,4:7); % for mediation analyses only

fdir = fullfile(SDL.out,SDL.data_type{1}(1:2)); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149};
if strcmp(SDL.data_type{1}(4:end),'Age')
    idx11 = strcmp(T.Group,'PTSD') & (T.Age<nanmedian(T.Age));  % group 11, PTSD & Age < 33 (median)
    idx12 = strcmp(T.Group,'PTSD') & (T.Age>=nanmedian(T.Age)); % group 12, PTSD & Age >= 33 (median)
    idx21 = strcmp(T.Group,'CONT') & (T.Age<nanmedian(T.Age));  % group 21, CONT & Age < 33 (median)
    idx22 = strcmp(T.Group,'CONT') & (T.Age>=nanmedian(T.Age)); % group 22, CONT & Age >= 33 (median)
elseif strcmp(SDL.data_type{1}(4:end),'Gender')
    idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'M');  % group 11, PTSD & Male
    idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'F');  % group 12, PTSD & Female
    idx21 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'M');  % group 21, CONT & Male
    idx22 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'F');  % group 22, CONT & Female
elseif strcmp(SDL.data_type{1}(4:end),'Dep')
%     try
%         idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'1');  % group 11, PTSD & depression
%         idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'0');  % group 12, PTSD & no-depression
%         idx21 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'1');  % group 21, CONT & depression
%         idx22 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'0');  % group 22, CONT & no-depression
%     catch
        idx11 = strcmp(T.Group,'PTSD') & (T.Dep==1);  % group 11, PTSD & depression
        idx12 = strcmp(T.Group,'PTSD') & (T.Dep==0);  % group 12, PTSD & no-depression
        idx21 = strcmp(T.Group,'CONT') & (T.Dep==1);  % group 21, CONT & depression
        idx22 = strcmp(T.Group,'CONT') & (T.Dep==0);  % group 22, CONT & no-depression
%     end
else
end
if strcmp(SDL.ana_type{1},'corr')
    func = @corr;
elseif strcmp(SDL.ana_type{1},'partialcorr')
    func = @partialcorr;
else
end
S0 = SDL_SCA4(R,idx11,idx12,idx21,idx22,func,0);     % SC of original PTSD & CONT data
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    X11(i-1) = SC0(i).mean11; X12(i-1) = SC0(i).mean12; X21(i-1) = SC0(i).mean21; X22(i-1) = SC0(i).mean22;
    
    p11(i-1) = SC1(i).p11; p12(i-1) = SC1(i).p12; p21(i-1) = SC1(i).p21; p22(i-1) = SC1(i).p22;
    
    X11R(i-1) = mean(SC1(i).mean11); X12R(i-1) = mean(SC1(i).mean12); X21R(i-1) = mean(SC1(i).mean21); X22R(i-1) = mean(SC1(i).mean22);
    X11RA(i-1,:) = SC1(i).mean11; X12RA(i-1,:) = SC1(i).mean12; X21RA(i-1,:) = SC1(i).mean21; X22RA(i-1,:) = SC1(i).mean22;
    
    CI11L(i-1) = SC1(i).CI11(1); CI11R(i-1) = SC1(i).CI11(2); CI12L(i-1) = SC1(i).CI12(1); CI12R(i-1) = SC1(i).CI12(2);
    CI21L(i-1) = SC1(i).CI21(1); CI21R(i-1) = SC1(i).CI21(2); CI22L(i-1) = SC1(i).CI22(1); CI22R(i-1) = SC1(i).CI22(2);

    diff_11v21L(i-1) = SC1(i).diffCI_11v21(1); diff_11v21R(i-1) = SC1(i).diffCI_11v21(2);
    diff_12v22L(i-1) = SC1(i).diffCI_12v22(1); diff_12v22R(i-1) = SC1(i).diffCI_12v22(2);
    diff_11v12L(i-1) = SC1(i).diffCI_11v12(1); diff_11v12R(i-1) = SC1(i).diffCI_11v12(2);
    diff_21v22L(i-1) = SC1(i).diffCI_21v22(1); diff_21v22R(i-1) = SC1(i).diffCI_21v22(2);
    diff_main1L(i-1)  = SC1(i).diffCI_main1(1); diff_main1R(i-1) = SC1(i).diffCI_main1(2);
    diff_main2L(i-1)  = SC1(i).diffCI_main2(1); diff_main2R(i-1) = SC1(i).diffCI_main2(2);
    diff_interL(i-1)  = SC1(i).diffCI_inter(1); diff_interR(i-1) = SC1(i).diffCI_inter(2);
    
    
    p_diff_11v21(i-1) = SC1(i).p_diff_11v21; 
    p_diff_12v22(i-1) = SC1(i).p_diff_12v22;
    p_diff_11v12(i-1) = SC1(i).p_diff_11v12;
    p_diff_21v22(i-1) = SC1(i).p_diff_21v22;
    p_diff_main1(i-1) = SC1(i).p_diff_main1;
    p_diff_main2(i-1) = SC1(i).p_diff_main2;
    p_diff_inter(i-1) = SC1(i).p_diff_inter;
end

figure;
a = [1:4]';
b = sum([(X11-X11R)',...
    (X12-X12R)',...
    (X21-X21R)',...
    (X22-X22R)'])';
b= bar(0,b);
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',16)
% xlabel('PTSD Severity','FontName', 'Arial','FontSize',14); 
ylabel('AUC','FontName', 'Arial','FontSize',16);
xticklabels({' '});
box off
b(1).FaceColor = 'r';
b(2).FaceColor = 'm';
b(3).FaceColor = 'b';
b(4).FaceColor = 'g';


figure;
% txt1 = 'Dep'; txt2 = 'NoDep';
% subplot(1,2,2);
X = 2:size(S0.Y11,1); 
% x2 = [X, fliplr(X)]; inBetween = [diff_interL, fliplr(diff_interR)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11-X12-X21+X22),'r-',X,(X11R-X12R-X21R+X22R),'b-','LineWidth',2); % interaction
plot(X,(X11-X11R),'r',X,(X12-X12R),'m',X,(X21-X21R),'b',X,(X22-X22R),'g','LineWidth',2); % four groups
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
set(gca, 'FontName', 'Arial'); set(gca,'FontSize',14)
xlabel('Top-N Regions','FontName', 'Arial','FontSize',16); 
set(gca,'XTick',0:10:150); 
% set(gca,'XTick',0:1:20); 
ylabel('SA-Based Structural Covariance','FontName', 'Arial','FontSize',16);%legend('95% CI','Areas with Reduction','Random Areas'); 
title('Real-Random','FontName', 'Arial','FontSize',16);
% legend('95% CI','Real Network','Random Network','FontName', 'Arial','FontSize',16);
legend('PTSD & Male','PTSD & Female','Healthy Male','Healthy Female','FontName', 'Arial','FontSize',16);
box off




savefig(fullfile(SDL.out,'Figure 5.fig'));


%% End
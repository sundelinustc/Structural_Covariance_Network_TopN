function Table6(SDL)

% calculate values for Table 6, which is to show the modification of Age on PTSD-related between group differences.
% (1) extract both global p-values and individual p-values from both CT-
% and SA-based outputs, and across three types of networks, i.e. atrophy,
% inflation and stable networks
% (2) FDR correction across 2 types of data (CT vs. SA) and 3 types of
% networks (atrophy, inflation vs. stable) and 8 between-group comparisons

%% Parameters
SDL.data_type = {'CT_Age10';'SA_Age10'}; % 2 types of data
SDL.net_type  = {'Reduction';'Increase';'Nodiff'}; % 3 types of networks

for idata = 1:size(SDL.data_type,1) % per data type
    for inet = 1:size(SDL.net_type,1) % per network type
        
        %% Loading information
        % Residuals
        SDL.out = fullfile(SDL.path,['Outputs_ComBat_',SDL.net_type{inet}]);
        fdir = fullfile(SDL.out,SDL.data_type{idata});
        fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{idata},'_corr_.mat']);
        load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
        R = T{:,2:149}; % residules
        
        idx11 = strcmp(T.Group,'PTSD') & (T.Age<10);             % group 11, PTSD & Age < 10
        idx12 = strcmp(T.Group,'PTSD') & (T.Age>=10 & T.Age<15); % group 12, PTSD & 10<=Age<15
        idx13 = strcmp(T.Group,'PTSD') & (T.Age>=15 & T.Age<20); % group 13, PTSD & 15<=Age<20
        idx14 = strcmp(T.Group,'PTSD') & (T.Age>=20 & T.Age<30); % group 14, PTSD & 20<=Age<30
        idx15 = strcmp(T.Group,'PTSD') & (T.Age>=30 & T.Age<40); % group 15, PTSD & 30<=Age<40
        idx16 = strcmp(T.Group,'PTSD') & (T.Age>=40 & T.Age<50); % group 16, PTSD & 40<=Age<50
        idx17 = strcmp(T.Group,'PTSD') & (T.Age>=50 & T.Age<60); % group 17, PTSD & 50<=Age<60
        idx18 = strcmp(T.Group,'PTSD') & (T.Age>=60);            % group 18, PTSD & 60<=Age
        
        idx21 = strcmp(T.Group,'CONT') & (T.Age<10);             % group 21, CONT & Age < 10
        idx22 = strcmp(T.Group,'CONT') & (T.Age>=10 & T.Age<15); % group 22, CONT & 10<=Age<15
        idx23 = strcmp(T.Group,'CONT') & (T.Age>=15 & T.Age<20); % group 23, CONT & 15<=Age<20
        idx24 = strcmp(T.Group,'CONT') & (T.Age>=20 & T.Age<30); % group 24, CONT & 20<=Age<30
        idx25 = strcmp(T.Group,'CONT') & (T.Age>=30 & T.Age<40); % group 25, CONT & 30<=Age<40
        idx26 = strcmp(T.Group,'CONT') & (T.Age>=40 & T.Age<50); % group 26, CONT & 40<=Age<50
        idx27 = strcmp(T.Group,'CONT') & (T.Age>=50 & T.Age<60); % group 27, CONT & 50<=Age<60
        idx28 = strcmp(T.Group,'CONT') & (T.Age>=60);            % group 28, CONT & 60<=Age
        
        % Effect size
        fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{idata},'_corr_.mat']);
        load(fn,'kk','d'); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);
        % T    --- table containing raw CT values and sbj info
        % R    --- residules of CT after lme regression
        % idx1 --- row index of PTSD in T & R
        % idx2 --- row index of CONT in T & R
        % d    --- Cohen's d for PTSD vs CONT
        % kk   --- the list of ranked area No. for PTSD vs CONT based on d, i.e. kk(1) showes the area with the largest CT reduction
        
        % Region's label & XYZ coordinates
        fdir = fullfile(SDL.path,'Original');
        load(fullfile(fdir,'RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations
        tbl = cell2table(D,'VariableNames',{'Area','X','Y','Z'});
        
        
        %% Structural Covariance Matrix -- Original
        S0 = SDL_SCA10(R,idx11,idx12,idx13,idx14,idx15,idx16,idx17,idx18,idx21,idx22,idx23,idx24,idx25,idx26,idx27,idx28,@corr,0);     % SC of original PTSD & CONT data
        % SC of original PTSD & CONT data
        % S0.Y1 --- R-to-Z transformation of corr/partialcorr matrix of group1
        % S0.Y2 --- R-to-Z transformation of corr/partialcorr matrix of group2
        
        %% Plot SC values of Top-N areas with changes in PTSD-CONT
        fdir = fullfile(SDL.out,SDL.data_type{idata});
        fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{idata},'_corr.mat']);
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
            
            
            
            
            X11RA(i-1,:) = SC1(i).mean11;
            X12RA(i-1,:) = SC1(i).mean12;
            X13RA(i-1,:) = SC1(i).mean13;
            X14RA(i-1,:) = SC1(i).mean14;
            X15RA(i-1,:) = SC1(i).mean15;
            X16RA(i-1,:) = SC1(i).mean16;
            X17RA(i-1,:) = SC1(i).mean17;
            X18RA(i-1,:) = SC1(i).mean18;
            
            X21RA(i-1,:) = SC1(i).mean21;
            X22RA(i-1,:) = SC1(i).mean22;
            X23RA(i-1,:) = SC1(i).mean23;
            X24RA(i-1,:) = SC1(i).mean24;
            X25RA(i-1,:) = SC1(i).mean25;
            X26RA(i-1,:) = SC1(i).mean26;
            X27RA(i-1,:) = SC1(i).mean27;
            X28RA(i-1,:) = SC1(i).mean28;
            
            %     CI11L(i-1) = SC1(i).CI11(1); CI11R(i-1) = SC1(i).CI11(2); CI12L(i-1) = SC1(i).CI12(1); CI12R(i-1) = SC1(i).CI12(2);
            %     CI21L(i-1) = SC1(i).CI21(1); CI21R(i-1) = SC1(i).CI21(2); CI22L(i-1) = SC1(i).CI22(1); CI22R(i-1) = SC1(i).CI22(2);
            %
            %     diff_11v21L(i-1) = SC1(i).diffCI_11v21(1); diff_11v21R(i-1) = SC1(i).diffCI_11v21(2);
            %     diff_12v22L(i-1) = SC1(i).diffCI_12v22(1); diff_12v22R(i-1) = SC1(i).diffCI_12v22(2);
            %     diff_11v12L(i-1) = SC1(i).diffCI_11v12(1); diff_11v12R(i-1) = SC1(i).diffCI_11v12(2);
            %     diff_21v22L(i-1) = SC1(i).diffCI_21v22(1); diff_21v22R(i-1) = SC1(i).diffCI_21v22(2);
            %     diff_main1L(i-1)  = SC1(i).diffCI_main1(1); diff_main1R(i-1) = SC1(i).diffCI_main1(2);
            %     diff_main2L(i-1)  = SC1(i).diffCI_main2(1); diff_main2R(i-1) = SC1(i).diffCI_main2(2);
            %     diff_interL(i-1)  = SC1(i).diffCI_inter(1); diff_interR(i-1) = SC1(i).diffCI_inter(2);
            %
            %
            p_diff_11v21(i-1) = SC1(i).p_diff_11v21;
            p_diff_12v22(i-1) = SC1(i).p_diff_12v22;
            p_diff_13v23(i-1) = SC1(i).p_diff_13v23;
            p_diff_14v24(i-1) = SC1(i).p_diff_14v24;
            p_diff_15v25(i-1) = SC1(i).p_diff_15v25;
            p_diff_16v26(i-1) = SC1(i).p_diff_16v26;
            p_diff_17v27(i-1) = SC1(i).p_diff_17v27;
            p_diff_18v28(i-1) = SC1(i).p_diff_18v28;
            %     p_diff_11v12(i-1) = SC1(i).p_diff_11v12;
            %     p_diff_21v22(i-1) = SC1(i).p_diff_21v22;
            %     p_diff_main1(i-1) = SC1(i).p_diff_main1;
            %     p_diff_main2(i-1) = SC1(i).p_diff_main2;
            %     p_diff_inter(i-1) = SC1(i).p_diff_inter;
        end
        
        
        
        %% individual values
        X11R=mean(X11RA,2)'; X21R=mean(X21RA,2)';
        X12R=mean(X12RA,2)'; X22R=mean(X22RA,2)';
        X13R=mean(X13RA,2)'; X23R=mean(X23RA,2)';
        X14R=mean(X14RA,2)'; X24R=mean(X24RA,2)';
        X15R=mean(X15RA,2)'; X25R=mean(X25RA,2)';
        X16R=mean(X16RA,2)'; X26R=mean(X26RA,2)';
        X17R=mean(X17RA,2)'; X27R=mean(X27RA,2)';
        X18R=mean(X18RA,2)'; X28R=mean(X28RA,2)';
        
        % mean of actual - random, 3(nets)x2(data)x8(between-groups)x147(comparisons)
        m_individual_diff(inet,idata,1,1,:) = X11-X11R; m_individual_diff(inet,idata,2,1,:) = X21-X21R;
        m_individual_diff(inet,idata,1,2,:) = X12-X12R; m_individual_diff(inet,idata,2,2,:) = X22-X22R;
        m_individual_diff(inet,idata,1,3,:) = X13-X13R; m_individual_diff(inet,idata,2,3,:) = X23-X23R;
        m_individual_diff(inet,idata,1,4,:) = X14-X14R; m_individual_diff(inet,idata,2,4,:) = X24-X24R;
        m_individual_diff(inet,idata,1,5,:) = X15-X15R; m_individual_diff(inet,idata,2,5,:) = X25-X25R;
        m_individual_diff(inet,idata,1,6,:) = X16-X16R; m_individual_diff(inet,idata,2,6,:) = X26-X26R;
        m_individual_diff(inet,idata,1,7,:) = X17-X17R; m_individual_diff(inet,idata,2,7,:) = X27-X27R;
        m_individual_diff(inet,idata,1,8,:) = X18-X18R; m_individual_diff(inet,idata,2,8,:) = X28-X28R;
        
        % p-values, random, 3(nets)x2(data)x8(between-group contrasts)x147(comparisons)
        p_individual_diff(inet,idata,1,:) = p_diff_11v21;
        p_individual_diff(inet,idata,2,:) = p_diff_12v22;
        p_individual_diff(inet,idata,3,:) = p_diff_13v23;
        p_individual_diff(inet,idata,4,:) = p_diff_14v24;
        p_individual_diff(inet,idata,5,:) = p_diff_15v25;
        p_individual_diff(inet,idata,6,:) = p_diff_16v26;
        p_individual_diff(inet,idata,7,:) = p_diff_17v27;
        p_individual_diff(inet,idata,8,:) = p_diff_18v28;
        
        
        % sign, random, 3(nets)x2(data)x10(between-group contrasts)x147(comparisons)
        sign_individual_diff(inet,idata,1,:) = sign((X11-X11R)-(X21-X21R));
        sign_individual_diff(inet,idata,2,:) = sign((X12-X12R)-(X22-X22R));
        sign_individual_diff(inet,idata,3,:) = sign((X13-X13R)-(X23-X23R));
        sign_individual_diff(inet,idata,4,:) = sign((X14-X14R)-(X24-X24R));
        sign_individual_diff(inet,idata,5,:) = sign((X15-X15R)-(X25-X25R));
        sign_individual_diff(inet,idata,6,:) = sign((X16-X16R)-(X26-X26R));
        sign_individual_diff(inet,idata,7,:) = sign((X17-X17R)-(X27-X27R));
        sign_individual_diff(inet,idata,8,:) = sign((X18-X18R)-(X28-X28R));
        
        
        %%  global values
        % mean of actual - random, 1x1
        m_global(inet,idata,1,1) = mean(X11-X11R); m_global(inet,idata,2,1) = mean(X21-X21R);
        m_global(inet,idata,1,2) = mean(X12-X12R); m_global(inet,idata,2,2) = mean(X22-X22R);
        m_global(inet,idata,1,3) = mean(X13-X13R); m_global(inet,idata,2,3) = mean(X23-X23R);
        m_global(inet,idata,1,4) = mean(X14-X14R); m_global(inet,idata,2,4) = mean(X24-X24R);
        m_global(inet,idata,1,5) = mean(X15-X15R); m_global(inet,idata,2,5) = mean(X25-X25R);
        m_global(inet,idata,1,6) = mean(X16-X16R); m_global(inet,idata,2,6) = mean(X26-X26R);
        m_global(inet,idata,1,7) = mean(X17-X17R); m_global(inet,idata,2,7) = mean(X27-X27R);
        m_global(inet,idata,1,8) = mean(X18-X18R); m_global(inet,idata,2,8) = mean(X28-X28R);
       
        
        
        % 95% CI (in fact 2.5~97.5% percentile), 1x2
        a = prctile(mean(X11'-X11RA),[2.5 97.5],2); err_global_L(inet,idata,1,1) = a(1); err_global_H(inet,idata,1,1) = a(2);
        a = prctile(mean(X21'-X21RA),[2.5 97.5],2); err_global_L(inet,idata,2,1) = a(1); err_global_H(inet,idata,2,1) = a(2);
        a = prctile(mean(X12'-X12RA),[2.5 97.5],2); err_global_L(inet,idata,1,2) = a(1); err_global_H(inet,idata,1,2) = a(2);
        a = prctile(mean(X22'-X22RA),[2.5 97.5],2); err_global_L(inet,idata,2,2) = a(1); err_global_H(inet,idata,2,2) = a(2);
        a = prctile(mean(X13'-X13RA),[2.5 97.5],2); err_global_L(inet,idata,1,3) = a(1); err_global_H(inet,idata,1,3) = a(2);
        a = prctile(mean(X23'-X23RA),[2.5 97.5],2); err_global_L(inet,idata,2,3) = a(1); err_global_H(inet,idata,2,3) = a(2);
        a = prctile(mean(X14'-X14RA),[2.5 97.5],2); err_global_L(inet,idata,1,4) = a(1); err_global_H(inet,idata,1,4) = a(2);
        a = prctile(mean(X24'-X24RA),[2.5 97.5],2); err_global_L(inet,idata,2,4) = a(1); err_global_H(inet,idata,2,4) = a(2);
        a = prctile(mean(X15'-X15RA),[2.5 97.5],2); err_global_L(inet,idata,1,5) = a(1); err_global_H(inet,idata,1,5) = a(2);
        a = prctile(mean(X25'-X25RA),[2.5 97.5],2); err_global_L(inet,idata,2,5) = a(1); err_global_H(inet,idata,2,5) = a(2);
        a = prctile(mean(X16'-X16RA),[2.5 97.5],2); err_global_L(inet,idata,1,6) = a(1); err_global_H(inet,idata,1,6) = a(2);
        a = prctile(mean(X26'-X26RA),[2.5 97.5],2); err_global_L(inet,idata,2,6) = a(1); err_global_H(inet,idata,2,6) = a(2);
        a = prctile(mean(X17'-X17RA),[2.5 97.5],2); err_global_L(inet,idata,1,7) = a(1); err_global_H(inet,idata,1,7) = a(2);
        a = prctile(mean(X27'-X27RA),[2.5 97.5],2); err_global_L(inet,idata,2,7) = a(1); err_global_H(inet,idata,2,7) = a(2);
        a = prctile(mean(X18'-X18RA),[2.5 97.5],2); err_global_L(inet,idata,1,8) = a(1); err_global_H(inet,idata,1,8) = a(2);
        a = prctile(mean(X28'-X28RA),[2.5 97.5],2); err_global_L(inet,idata,2,8) = a(1); err_global_H(inet,idata,2,8) = a(2);
        
        
        % p values, 1x1
        if mean(X11-X21)>0, p_global_diff_11v21(inet,idata)=SDL_p_permutation(mean(X11-X21),mean(X11RA-X21RA,1)); else,p_global_diff_11v21(inet,idata)=SDL_p_permutation(mean(X21-X11),mean(X21RA-X11RA,1));end
        if mean(X12-X22)>0, p_global_diff_12v22(inet,idata)=SDL_p_permutation(mean(X12-X22),mean(X12RA-X22RA,1)); else,p_global_diff_12v22(inet,idata)=SDL_p_permutation(mean(X22-X12),mean(X22RA-X12RA,1));end
        if mean(X13-X23)>0, p_global_diff_13v23(inet,idata)=SDL_p_permutation(mean(X13-X23),mean(X13RA-X23RA,1)); else,p_global_diff_13v23(inet,idata)=SDL_p_permutation(mean(X23-X13),mean(X23RA-X13RA,1));end
        if mean(X14-X24)>0, p_global_diff_14v24(inet,idata)=SDL_p_permutation(mean(X14-X24),mean(X14RA-X24RA,1)); else,p_global_diff_14v24(inet,idata)=SDL_p_permutation(mean(X24-X14),mean(X24RA-X14RA,1));end
        if mean(X15-X25)>0, p_global_diff_15v25(inet,idata)=SDL_p_permutation(mean(X15-X25),mean(X15RA-X25RA,1)); else,p_global_diff_15v25(inet,idata)=SDL_p_permutation(mean(X25-X15),mean(X25RA-X15RA,1));end
        if mean(X16-X26)>0, p_global_diff_16v26(inet,idata)=SDL_p_permutation(mean(X16-X26),mean(X16RA-X26RA,1)); else,p_global_diff_16v26(inet,idata)=SDL_p_permutation(mean(X26-X16),mean(X26RA-X16RA,1));end
        if mean(X17-X27)>0, p_global_diff_17v27(inet,idata)=SDL_p_permutation(mean(X17-X27),mean(X17RA-X27RA,1)); else,p_global_diff_17v27(inet,idata)=SDL_p_permutation(mean(X27-X17),mean(X27RA-X17RA,1));end
        if mean(X18-X28)>0, p_global_diff_18v28(inet,idata)=SDL_p_permutation(mean(X18-X28),mean(X18RA-X28RA,1)); else,p_global_diff_18v28(inet,idata)=SDL_p_permutation(mean(X28-X18),mean(X28RA-X18RA,1));end

        
        
        % sign
        sign_global_diff_11v21(inet,idata) = sign(mean((X11-X11R)-(X21-X21R)));
        sign_global_diff_12v22(inet,idata) = sign(mean((X12-X12R)-(X22-X22R)));
        sign_global_diff_13v23(inet,idata) = sign(mean((X13-X13R)-(X23-X23R)));
        sign_global_diff_14v24(inet,idata) = sign(mean((X14-X14R)-(X24-X24R)));
        sign_global_diff_15v25(inet,idata) = sign(mean((X15-X15R)-(X25-X25R)));
        sign_global_diff_16v26(inet,idata) = sign(mean((X16-X16R)-(X26-X26R)));
        sign_global_diff_17v27(inet,idata) = sign(mean((X17-X17R)-(X27-X27R)));
        sign_global_diff_18v28(inet,idata) = sign(mean((X18-X18R)-(X28-X28R)));

        
        
    end
end


%% plot
figure
for i1 = 1:3
    for i2 = 1:2
        subplot(3,2,(i1-1)*2+i2);
        
        y=[]; y = reshape(m_global(i1,i2,1,:),8,1); y = [y,reshape(m_global(i1,i2,2,:),8,1)];%y=147*y;
        errlow  = reshape(err_global_L(i1,i2,1,:),8,1); errlow  = [errlow, reshape(err_global_L(i1,i2,2,:),8,1)];%errlow=147*errlow;
        errhigh = reshape(err_global_H(i1,i2,1,:),8,1); errhigh = [errhigh,reshape(err_global_H(i1,i2,2,:),8,1)];%errhigh=147*errhigh;
%         model_error = (errhigh-errlow)/2;
        
        b = bar(y, 'grouped'); b(1).FaceColor = 'r'; b(2).FaceColor = 'b';
        ylabel('Adj. Mean SC');
        set(gca,'xticklabel',{'<10','10~15','15~20','20~30','30~40','40~50','50~60','>=60'}); 
        set(gca,'box','off');
        xtickangle(45);
        hold on %For MATLAB 2019b or later releases
        
        nbars = size(y, 2);% Calculate the number of bars in each group
        x = []; % Get the x coordinate of the bars
        for i = 1:nbars, x = [x ; b(i).XEndPoints]; end
        errorbar(x',y,y-errlow,errhigh-y,'k','linestyle','none')'; % Plot the errorbars
%         errorbar(x',y,model_error,'k','linestyle','none')'; % Plot the errorbars
        hold off
        
        
    end
end


%% Global p-values (FDR corrected across 2 data type, 3 network types and 6 comparisons)
p_global_diff_all(:,:,1) = p_global_diff_11v21; sign_global_diff_all(:,:,1) =  sign_global_diff_11v21;
p_global_diff_all(:,:,2) = p_global_diff_12v22; sign_global_diff_all(:,:,2) =  sign_global_diff_12v22;
p_global_diff_all(:,:,3) = p_global_diff_13v23; sign_global_diff_all(:,:,3) =  sign_global_diff_13v23;
p_global_diff_all(:,:,4) = p_global_diff_14v24; sign_global_diff_all(:,:,4) =  sign_global_diff_14v24;
p_global_diff_all(:,:,5) = p_global_diff_15v25; sign_global_diff_all(:,:,5) =  sign_global_diff_15v25;
p_global_diff_all(:,:,6) = p_global_diff_16v26; sign_global_diff_all(:,:,6) =  sign_global_diff_16v26;
p_global_diff_all(:,:,7) = p_global_diff_17v27; sign_global_diff_all(:,:,7) =  sign_global_diff_17v27;
p_global_diff_all(:,:,8) = p_global_diff_18v28; sign_global_diff_all(:,:,8) =  sign_global_diff_18v28;

[~, ~, ~, p_global_diff_all_adj] = fdr_bh(p_global_diff_all,0.05,'pdep','yes');
p_global_diff_all_adj(p_global_diff_all_adj<0.0001) = 0.0001; % to replace the values == 0
p_global_diff_all_adj = p_global_diff_all_adj.*sign_global_diff_all; % p-values with sign

reshape(p_global_diff_all_adj(1,1,:),1,8) 
reshape(p_global_diff_all_adj(1,2,:),1,8)

% for idata = 1:size(SDL.data_type,1) % per data type
%     for inet = 1:size(SDL.net_type,1) % per network type
%         p_global_diff(inet,idata,1:4,1:4) = 1; % initial values, NS
%         p_global_diff(inet,idata,1,2) = p_global_diff_all_adj(inet,idata,1); p_global_diff(inet,idata,2,1) = p_global_diff_all_adj(inet,idata,1);
%         p_global_diff(inet,idata,1,3) = p_global_diff_all_adj(inet,idata,2); p_global_diff(inet,idata,3,1) = p_global_diff_all_adj(inet,idata,2);
%         p_global_diff(inet,idata,1,4) = p_global_diff_all_adj(inet,idata,3); p_global_diff(inet,idata,4,1) = p_global_diff_all_adj(inet,idata,3);
%         p_global_diff(inet,idata,1,5) = p_global_diff_all_adj(inet,idata,4); p_global_diff(inet,idata,5,1) = p_global_diff_all_adj(inet,idata,4);
%         p_global_diff(inet,idata,2,3) = p_global_diff_all_adj(inet,idata,5); p_global_diff(inet,idata,3,2) = p_global_diff_all_adj(inet,idata,5);
%         p_global_diff(inet,idata,2,4) = p_global_diff_all_adj(inet,idata,6); p_global_diff(inet,idata,4,2) = p_global_diff_all_adj(inet,idata,6);
%         p_global_diff(inet,idata,2,5) = p_global_diff_all_adj(inet,idata,7); p_global_diff(inet,idata,5,2) = p_global_diff_all_adj(inet,idata,7);
%         p_global_diff(inet,idata,3,4) = p_global_diff_all_adj(inet,idata,8); p_global_diff(inet,idata,4,3) = p_global_diff_all_adj(inet,idata,8);
%         p_global_diff(inet,idata,3,5) = p_global_diff_all_adj(inet,idata,9); p_global_diff(inet,idata,5,3) = p_global_diff_all_adj(inet,idata,9);
%         p_global_diff(inet,idata,4,5) =p_global_diff_all_adj(inet,idata,10); p_global_diff(inet,idata,5,4) = p_global_diff_all_adj(inet,idata,10);
%     end
% end


%% Individual p-values (FDR corrected across 2 data type, 3 network types, 6 contrasts and 147 tests)

[~, ~, ~, p_individual_diff_adj] = fdr_bh(p_individual_diff,0.05,'pdep','yes');
% a, the matrix to show both significant tests (P<0.05) and the sign of the tests (e.g. PTSD>random, or PTSD>nonPTSD)
% J, index of the significant tests
% J+1 because tests 1:147 corresponding to top-2 to -148 regions
% b is 1x1x107 for example, and b(:) is 107x1x1, consistent with J

p_individual_diff_adj(p_individual_diff_adj<0.0001) = 0.0001; % to replace the values == 0
fprintf('\n===Areas,individual_diff,FDR corrected\n');
a = (p_individual_diff_adj<=0.05).* sign_individual_diff;
for idata = 1:size(SDL.data_type,1)
    for inet = 1:size(SDL.net_type,1)
        for icon = 1:10
            fprintf('data=%s, net=%s, con=%d, areas: ',SDL.data_type{idata},SDL.net_type{inet},icon);
            J = find(a(inet,idata,icon,:));
            b = a(inet,idata,icon,J);
            fprintf('%d ',(J+1).*b(:));fprintf('\n');
        end
    end
end


% END
end
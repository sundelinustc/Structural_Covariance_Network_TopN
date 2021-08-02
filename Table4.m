function Table4(SDL)

% calculate values for Table 4, which is to show the PTSD x depression interaction.
% (1) extract both global p-values and individual p-values from both CT-
% and SA-based outputs, and across three types of networks, i.e. atrophy,
% inflation and stable networks
% (2) FDR correction across 2 types of data (CT vs. SA) and 3 types of
% networks (atrophy, inflation vs. stable) and 6 between-group comparisons

%% Parameters
SDL.data_type = {'CT_Dep';'SA_Dep'}; % 2 types of data
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
        
        idx11 = strcmp(T.Group,'PTSD') & (T.Dep==1);  % group 11, PTSD & depression
        idx12 = strcmp(T.Group,'PTSD') & (T.Dep==0);  % group 12, PTSD & no-depression
        idx21 = strcmp(T.Group,'CONT') & (T.Dep==1);  % group 21, CONT & depression
        idx22 = strcmp(T.Group,'CONT') & (T.Dep==0);  % group 22, CONT & no-depression
        
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
        S0 = SDL_SCA4(R,idx11,idx12,idx21,idx22,@corr,0);     % SC of original PTSD & CONT data
        % S0.Y1 --- R-to-Z transformation of corr/partialcorr matrix of group1
        % S0.Y2 --- R-to-Z transformation of corr/partialcorr matrix of group2
        
        %% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT
        fdir = fullfile(SDL.out,SDL.data_type{idata});
        fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{idata},'_corr.mat']);
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
        
        %% individual values
        % mean of actual - random, 3(nets)x2(data)x4(groups)x147(comparisons)
        m_individual_diff(inet,idata,1,:) = X11-X11R; 
        m_individual_diff(inet,idata,2,:) = X21-X21R;
        m_individual_diff(inet,idata,3,:) = X12-X12R;
        m_individual_diff(inet,idata,4,:) = X22-X22R;
        
        % p-values, random, 3(nets)x2(data)x4(between-group contrasts)x147(comparisons)
        p_individual_diff(inet,idata,1,:) = p_diff_11v21; 
        p_individual_diff(inet,idata,2,:) = p_diff_12v22; 
        p_individual_diff(inet,idata,3,:) = p_diff_11v12; 
        p_individual_diff(inet,idata,4,:) = p_diff_21v22; 
        p_individual_inter(inet,idata,:)  = p_diff_inter; % interaction, 3x2x147

        
        % sign, random, 3(nets)x2(data)x10(between-group contrasts)x147(comparisons)
        sign_individual_diff(inet,idata,1,:) = sign((X11-X11R)-(X21-X21R)); 
        sign_individual_diff(inet,idata,2,:) = sign((X12-X12R)-(X22-X22R));
        sign_individual_diff(inet,idata,3,:) = sign((X11-X11R)-(X12-X12R));
        sign_individual_diff(inet,idata,4,:) = sign((X21-X21R)-(X22-X22R));
        
        
        %%  global values
        % mean of actual - random, 1x1
        m_global(inet,idata,1,1) = mean(X11-X11R);
        m_global(inet,idata,2,1) = mean(X21-X21R);
        m_global(inet,idata,1,2) = mean(X12-X12R);
        m_global(inet,idata,2,2) = mean(X22-X22R);

        
        % 95% CI (in fact 2.5~97.5% percentile), 1x2
        a = prctile(mean(X11'-X11RA),[2.5 97.5],2); err_global_L(inet,idata,1,1) = a(1); err_global_H(inet,idata,1,1) = a(2);
        a = prctile(mean(X21'-X21RA),[2.5 97.5],2); err_global_L(inet,idata,2,1) = a(1); err_global_H(inet,idata,2,1) = a(2);
        a = prctile(mean(X12'-X12RA),[2.5 97.5],2); err_global_L(inet,idata,1,2) = a(1); err_global_H(inet,idata,1,2) = a(2);
        a = prctile(mean(X22'-X22RA),[2.5 97.5],2); err_global_L(inet,idata,2,2) = a(1); err_global_H(inet,idata,2,2) = a(2);
                
        
        % p values, 1x1
        if mean(X11-X21)>0, p_global_diff_11v21(inet,idata)=SDL_p_permutation(mean(X11-X21),mean(X11RA-X21RA,1)); else,p_global_diff_11v21(inet,idata)=SDL_p_permutation(mean(X21-X11),mean(X21RA-X11RA,1));end
        if mean(X12-X22)>0, p_global_diff_12v22(inet,idata)=SDL_p_permutation(mean(X12-X22),mean(X12RA-X22RA,1)); else,p_global_diff_12v22(inet,idata)=SDL_p_permutation(mean(X22-X12),mean(X22RA-X12RA,1));end
        if mean(X11-X12)>0, p_global_diff_11v12(inet,idata)=SDL_p_permutation(mean(X11-X12),mean(X11RA-X12RA,1)); else,p_global_diff_11v12(inet,idata)=SDL_p_permutation(mean(X12-X11),mean(X12RA-X11RA,1));end
        if mean(X21-X22)>0, p_global_diff_21v22(inet,idata)=SDL_p_permutation(mean(X21-X22),mean(X21RA-X22RA,1)); else,p_global_diff_21v22(inet,idata)=SDL_p_permutation(mean(X22-X21),mean(X22RA-X21RA,1));end
        if mean(X12-X21)>0, p_global_diff_12v21(inet,idata)=SDL_p_permutation(mean(X12-X21),mean(X12RA-X21RA,1)); else,p_global_diff_12v21(inet,idata)=SDL_p_permutation(mean(X21-X12),mean(X21RA-X12RA,1));end
        if mean(X21-X22)>0, p_global_diff_11v22(inet,idata)=SDL_p_permutation(mean(X11-X22),mean(X11RA-X22RA,1)); else,p_global_diff_11v22(inet,idata)=SDL_p_permutation(mean(X22-X11),mean(X22RA-X11RA,1));end

        if mean((X11-X21)-(X12-X22))>0
            p_global_inter(inet,idata)=SDL_p_permutation(mean((X11-X21)-(X12-X22)),mean((X11RA-X21RA)-(X12RA-X22RA),1)); 
        else
            p_global_inter(inet,idata)=SDL_p_permutation(mean((X12-X22)-(X11-X21)),mean((X12RA-X22RA)-(X11RA-X21RA),1)); 
        end
             
        % sign
        sign_global_diff_11v21(inet,idata) = sign(mean((X11-X11R)-(X21-X21R)));
        sign_global_diff_12v22(inet,idata) = sign(mean((X12-X12R)-(X22-X22R)));
        sign_global_diff_11v12(inet,idata) = sign(mean((X11-X11R)-(X12-X12R)));
        sign_global_diff_21v22(inet,idata) = sign(mean((X21-X21R)-(X22-X22R)));
        sign_global_diff_12v21(inet,idata) = sign(mean((X12-X12R)-(X21-X21R)));
        sign_global_diff_11v22(inet,idata) = sign(mean((X11-X11R)-(X22-X22R)));
        
        
    end
end


%% plot
figure
for i1 = 1:3
    for i2 = 1:2
        subplot(3,2,(i1-1)*2+i2);
        x = 1:4; 
        y = reshape(m_global(i1,i2,:),1,4); 
        errlow  = reshape(err_global_L(i1,i2,:),1,4); 
        errhigh = reshape(err_global_H(i1,i2,:),1,4);
        bar(x,y); 
        ylabel('Adj. Mean SC');
        set(gca,'xticklabel',{'PTSD&Dep','Dep','PTSD','Control'}); 
        set(gca,'box','off'); 
        hold on; 
        er = errorbar(x,y,y-errlow,errhigh-y); 
        er.Color = [0 0 0]; 
        er.LineStyle = 'none'; hold off;
    end
end


%% Global p-values (FDR corrected across 2 data type, 3 network types and 6 comparisons)
[~, ~, ~, p_global_inter_adj] = fdr_bh(p_global_inter,0.05,'pdep','yes')


p_global_diff_all(:,:,1) = p_global_diff_11v21; sign_global_diff_all(:,:,1) =  sign_global_diff_11v21;
p_global_diff_all(:,:,2) = p_global_diff_12v22; sign_global_diff_all(:,:,2) =  sign_global_diff_12v22;
p_global_diff_all(:,:,3) = p_global_diff_11v12; sign_global_diff_all(:,:,3) =  sign_global_diff_11v12;
p_global_diff_all(:,:,4) = p_global_diff_21v22; sign_global_diff_all(:,:,4) =  sign_global_diff_21v22;
p_global_diff_all(:,:,5) = p_global_diff_12v21; sign_global_diff_all(:,:,5) =  sign_global_diff_12v21;
p_global_diff_all(:,:,6) = p_global_diff_11v22; sign_global_diff_all(:,:,6) =  sign_global_diff_11v22;

[~, ~, ~, p_global_diff_all_adj] = fdr_bh(p_global_diff_all,0.05,'pdep','yes');
p_global_diff_all_adj = p_global_diff_all_adj.*sign_global_diff_all; % p-values with sign

reshape(p_global_diff_all_adj(1,2,:),1,6) % SA-based atrophy networks
reshape(p_global_diff_all_adj(3,2,:),1,6) % SA-based stable networks
  
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
% corrected for 3 (types of networks) x 2 (types of data) x 6 (between-group contrasts) x 147 (individual size)
[~, ~, ~, p_individual_inter_adj] = fdr_bh(p_individual_inter,0.05,'pdep','yes');

[~, ~, ~, p_individual_diff_adj] = fdr_bh(p_individual_diff,0.05,'pdep','yes');
% a, the matrix to show both significant tests (P<0.05) and the sign of the tests (e.g. PTSD>random, or PTSD>nonPTSD)
% J, index of the significant tests
% J+1 because tests 1:147 corresponding to top-2 to -148 regions
% b is 1x1x107 for example, and b(:) is 107x1x1, consistent with J

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
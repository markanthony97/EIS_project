load('continuastasera.mat')
predict = (y>0.48); % best threshold
rows = cell(59,2);
predicted = cell(59,2);
adult_numeric_pca = cell(59,2);
adult_cate_names = adult_ohe_names(4:63);
h = zeros(59,4);
p = zeros(59,4);
numeric_norma=array2table(zscore(table2array(adult_numeric)),'VariableNames',adult_numeric.Properties.VariableNames);
[~,numeric_pca,~] = pca(table2array(adult_numeric));
for i = 5:63
    rows(i-4,1)=adult_ohe_names(i);
    truth(i-4,1)=adult_ohe_names(i);
    predicted(i-4,1) = adult_ohe_names(i);
    adult_numeric_pca(i-4,1) = adult_ohe_names(i);
    rows(i-4,2)={boolean(adult_for_nn_reduced(:,i))};
    truth(i-4,2)={result(boolean(adult_for_nn_reduced(:,i)))};
    predicted(i-4,2) = {predict(boolean(adult_for_nn_reduced(:,i)))};
    % adult_numeric_pca(i-4,2) = {numeric_pca(boolean(adult_for_nn_reduced(:,i)),:)};
    for j=1:4
        [h(i,j),p(i,j)] = ttest(adult_for_nn_reduced(:,j),y');
    end
end

% save plotthesetowmorrow.mat rows truth predicted  adult_numeric_cat
%save pcatest.mat adult_numeric_pca 
%clear
%% 
% load plotthesettomorrow.mat
% unique_values
% NOW TRY GETTING 
% CREATE 8 ERRORBAR GRAPHS, ONE FOR EACH CATEGORICAL FEATURE, CONTAINING 
% THE DISTRIBUTION IN TERMS OF BAR GRAPH OF EACH VALUE
% https://www.mathworks.com/help/releases/R2020b/matlab/ref/barh.html
% https://www.mathworks.com/help/releases/R2020b/matlab/ref/bar3h.html
clear
load plotthesetowmorrow.mat
load pcatest.mat
load grouped_clusters.mat
grouped_clusters.native_country(grouped_clusters.native_country~='United-States') = 'Foreign-Country';
variables = grouped_clusters(:,5:12).Properties.VariableNames;
variables{3} = "marital\_status";
variables{8} = "native\_country";
cmap = hsv;
for j=1:8
    f1 = figure(j);
    Ax(1) = axes(f1); 
    %set(Ax(1), 'XScale', 'log')
    set(Ax(1), 'ZScale', 'log')
    %set(Ax(1), 'ZTicks', 'log')
    unique_vals{j} = unique(table2array(grouped_clusters(:,j+4)));
    num_cats = length(unique_vals{j});
    xlabel('age');
    ylabel('hours per week');
    zlabel('capital gains')
    title(strcat("Groups of ",variables{j}))
    hold on
    randcolors = cmap(1:floor(255/num_cats):255,:);
    
    for i=1:num_cats
       cluster_index=find(rows(:,1)==unique_vals{j}(i));
       subtable = adult_numeric_cat{cluster_index,2};
       
       %subtable.net_yearly_capital = subtable.capital_gain - subtable.capital_loss;
       subtable.capital_gain(subtable.capital_gain==0)=1;
       %subtable.age = subtable.age;
       %subtable.hours_per_week = subtable.hours_per_week;
       subpredicted = predicted{cluster_index,2};
       subtruth = truth{cluster_index,2};
   subtable_TP = subtable(boolean(subtruth' .* subpredicted),:);
   subtable_TN = subtable(boolean(~subtruth' .* ~subpredicted),:);
   subtable_FP = subtable(boolean((subtruth==false)' .* (subpredicted==true)),:);
   subtable_FN = subtable(boolean((subtruth==true)' .* (subpredicted==false)),:);
       alpha = 0.9-(height(subtable))/height(grouped_clusters); 
       
       
       
       corralpha = alpha + 0.1;
       if corralpha>1
           corralpha=1.0;
       end
       if alpha<0.1
           alpha=alpha+0.15;
       end
%        alpha = 0.9;
%        corralpha= 1;
%        scatterTP(i)=scatter(subtable_TP(:,1),subtable_TP(:,2),15,'ks','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
%        scatterTN(i)=scatter(subtable_TN(:,1),subtable_TN(:,2),15,'kd','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
%        scatterFP(i)=scatter(subtable_FP(:,1),subtable_FP(:,2),35,'ko','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
%        scatterFN(i)=scatter(subtable_FN(:,1),subtable_FN(:,2),35,'k^','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       
       scatterTP(i)=scatter3(subtable_TP.age,subtable_TP.hours_per_week,subtable_TP.capital_gain,15,'ks','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       scatterTN(i)=scatter3(subtable_TN.age,subtable_TN.hours_per_week,subtable_TN.capital_gain,15,'kd','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       scatterFP(i)=scatter3(subtable_FP.age,subtable_FP.hours_per_week,subtable_FP.capital_gain,35,'ko','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       scatterFN(i)=scatter3(subtable_FN.age,subtable_FN.hours_per_week,subtable_FN.capital_gain,35,'k^','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
        
    end
    for i=1:num_cats
        legend_names(1+(i-1)*4:i*4) = unique_vals{j}(i);
        legend_names_2(1+(i-1)*4:i*4) = {'True positives','True negatives',
            'False positives','False negatives'};
    end
    legend_values = legend(Ax(1),scatterFP,legend_names(1:4:num_cats*4),'Location','northeastoutside');
    title(legend_values,'possible values')
    % view(Ax(1),[0.297068403908789 90]);
    view(Ax(1),[-102.339828249926 39.817358490566]);
    Ax(2) = copyobj(Ax(1),gcf);
    delete(get(Ax(2), 'Children'))
    hold on
    scatterEMPTY(1) = scatter3([],[],[],'ks','filled','Parent', Ax(2));
    scatterEMPTY(2) = scatter3([],[],[],'kd','filled','Parent', Ax(2));
    scatterEMPTY(3) = scatter3([],[],[],'ko','filled','Parent', Ax(2));
    scatterEMPTY(4) = scatter3([],[],[],'k^','filled','Parent', Ax(2));
    hold on
    set(Ax(2), 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    legend_predvstruth = legend(Ax(2), scatterEMPTY,legend_names_2(1:4),'Location','southeastoutside','Position',[0.838723769369846 0.143794162826421 0.0995607594925348 0.137480795108778]);
    title(legend_predvstruth,'comparison of values')
    grid
    hold off
    if j==3
        variables{3} = "marital_status";
    end
    if j==8
        variables{8} = "native_country";
    end
      %saveas(f1,strcat('scatterplots/scatterplot_cap_gain',variables{j},'.fig'))
      f1.Position=[1,20,1366,768];
      saveas(f1,strcat('scatterplots/scatterplot_cap_gain',variables{j},'.svg'))
end
%% 
clear
close all
load plotthesetowmorrow.mat
load pcatest.mat
load grouped_clusters.mat
grouped_clusters.native_country(grouped_clusters.native_country~='United-States') = 'Foreign-Country';
variables = grouped_clusters(:,5:12).Properties.VariableNames;
variables{3} = "marital\_status";
variables{8} = "native\_country";
cmap = hsv;

for j=1:8
    f1 = figure(j);
    Ax(1) = axes(f1); 
    %set(Ax(1), 'XScale', 'log')
    set(Ax(1), 'ZScale', 'log')
    %set(Ax(1), 'ZTicks', 'log')
    unique_vals{j} = unique(table2array(grouped_clusters(:,j+4)));
    num_cats = length(unique_vals{j});
    xlabel('age');
    ylabel('hours per week');
    zlabel('capital loss')
    title(strcat("Groups of ",variables{j}),'considering capital losses')
    hold on
    randcolors = cmap(1:floor(255/num_cats):255,:);
    
    for i=1:num_cats
       cluster_index=find(rows(:,1)==unique_vals{j}(i));
       subtable = adult_numeric_cat{cluster_index,2};
       
       %subtable.net_yearly_capital = subtable.capital_gain - subtable.capital_loss;
       subtable.capital_loss(subtable.capital_loss==0)=1;
       %subtable.age = subtable.age;
       %subtable.hours_per_week = subtable.hours_per_week;
       subpredicted = predicted{cluster_index,2};
       subtruth = truth{cluster_index,2};
   subtable_TP = subtable(boolean(subtruth' .* subpredicted),:);
   subtable_TN = subtable(boolean(~subtruth' .* ~subpredicted),:);
   subtable_FP = subtable(boolean((subtruth==false)' .* (subpredicted==true)),:);
   subtable_FN = subtable(boolean((subtruth==true)' .* (subpredicted==false)),:);
       alpha = 0.9-(height(subtable))/height(grouped_clusters); 
       
       
       
       corralpha = alpha + 0.1;
       if corralpha>1
           corralpha=1.0;
       end
       if alpha<0.1
           alpha=alpha+0.15;
       end

       scatterTP(i)=scatter3(subtable_TP.age,subtable_TP.hours_per_week,subtable_TP.capital_loss,15,'ks','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       scatterTN(i)=scatter3(subtable_TN.age,subtable_TN.hours_per_week,subtable_TN.capital_loss,15,'kd','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       scatterFP(i)=scatter3(subtable_FP.age,subtable_FP.hours_per_week,subtable_FP.capital_loss,35,'ko','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
       scatterFN(i)=scatter3(subtable_FN.age,subtable_FN.hours_per_week,subtable_FN.capital_loss,35,'k^','MarkerFaceColor',randcolors(i,:),'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',corralpha);
        
    end
    for i=1:num_cats
        legend_names(1+(i-1)*4:i*4) = unique_vals{j}(i);
        legend_names_2(1+(i-1)*4:i*4) = {'True positives','True negatives',
            'False positives','False negatives'};
    end
    legend_values = legend(Ax(1),scatterFP,legend_names(1:4:num_cats*4),'Location','northeastoutside');
    title(legend_values,'possible values')
    % view(Ax(1),[0.297068403908789 90]);
    view(Ax(1),[-102.339828249926 39.817358490566]);
    Ax(2) = copyobj(Ax(1),gcf);
    delete(get(Ax(2), 'Children'))
    hold on
    scatterEMPTY(1) = scatter3([],[],[],'ks','filled','Parent', Ax(2));
    scatterEMPTY(2) = scatter3([],[],[],'kd','filled','Parent', Ax(2));
    scatterEMPTY(3) = scatter3([],[],[],'ko','filled','Parent', Ax(2));
    scatterEMPTY(4) = scatter3([],[],[],'k^','filled','Parent', Ax(2));
    hold on
    set(Ax(2), 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    legend_predvstruth = legend(Ax(2), scatterEMPTY,legend_names_2(1:4),'Location','southeastoutside','Position',[0.838723769369846 0.143794162826421 0.0995607594925348 0.137480795108778]);
    title(legend_predvstruth,'comparison of values')
    grid
    hold off
    if j==3
        variables{3} = "marital_status";
    end
    if j==8
        variables{8} = "native_country";
    end
      %saveas(f1,strcat('scatterplots/scatterplot_cap_loss',variables{j},'.fig'))
      f1.Position=[1,20,1366,768];
      saveas(f1,strcat('scatterplots/scatterplot_cap_loss',variables{j},'.svg'))
end

% close all
%% numeric data statistics for each case
load plotthesetowmorrow.mat
% confusion_matrix = table(zeros(59,4),'RowNames',adult_numeric_cat(:,1));
for i=1:59
   subtable = adult_numeric_cat{i,2};
   subpredicted = predicted{i,2};
   subtruth = truth{i,2};
   subtable_TP = subtable(boolean(subtruth' .* subpredicted),:);
   subtable_TN = subtable(boolean(~subtruth' .* ~subpredicted),:);
   subtable_FP = subtable(boolean((subtruth==false)' .* (subpredicted==true)),:);
   subtable_FN = subtable(boolean((subtruth==true)' .* (subpredicted==false)),:);
   TP(i) = height(subtable_TP);
   FP(i) = height(subtable_FP);
   TN(i) = height(subtable_TN);
   FN(i) = height(subtable_FN);
   % capital gain statistics
   cap_gain_mean(i) = mean(subtable.capital_gain);
   cap_gain_std(i) = std(subtable.capital_gain);
   cap_gain_min(i) = min(subtable.capital_gain);
   cap_gain_max(i) = max(subtable.capital_gain);
   % capital loss statistics
   cap_loss_mean(i) = mean(subtable.capital_loss);
   cap_loss_std(i) = std(subtable.capital_loss);
   cap_loss_min(i) = min(subtable.capital_loss);
   cap_loss_max(i) = max(subtable.capital_loss);
   % age statistics
   age_mean(i) = mean(subtable.age);
   age_std(i) = std(subtable.age);
   age_min(i) = min(subtable.age);
   age_max(i) = max(subtable.age);
   % hours_per_week_statistics
   hpw_mean(i) = mean(subtable.hours_per_week);
   hpw_std(i) = std(subtable.hours_per_week);
   hpw_min(i) = min(subtable.hours_per_week);
   hpw_max(i) = min(subtable.hours_per_week);
end
cmat = table(TP',TN',FP',FN','VariableNames',{'TP','TN','FP','FN'},'RowNames',adult_numeric_cat(:,1));
cmat.accuracy = round(100*(cmat.TP+cmat.TN)./(cmat.TP+cmat.TN+cmat.FP+cmat.FN),2);
cmat.precision = round(100*cmat.TP./(cmat.TP+cmat.FP),2);
cmat.recall = round(100*cmat.TP./(cmat.TP+cmat.FN),2);

cmat_workclass = cmat(1:7,5:7);
cmat_education = cmat(8:23,5:7);
cmat_mar_stat = cmat(24:30,5:7);
cmat_occupation = cmat(31:44,5:7);
cmat_relationship = cmat(45:50,5:7);
cmat_race = cmat(51:55,5:7);
cmat_gender = cmat(56:57,5:7);
cmat_nc = cmat(58:59,5:7);

table2latex(cmat_workclass,'latextables/cmat_workclass')
table2latex(cmat_education,'latextables/cmat_education')
table2latex(cmat_mar_stat,'latextables/cmat_mar_stat')
table2latex(cmat_occupation,'latextables/cmat_occupation')
table2latex(cmat_relationship,'latextables/cmat_relationship')
table2latex(cmat_race,'latextables/cmat_race')
table2latex(cmat_gender,'latextables/cmat_gender')
table2latex(cmat_nc,'latextables/cmat_nc')


% round everything
cap_gain_mean = round(cap_gain_mean,2);
cap_gain_std = round(cap_gain_std,2);
cap_gain_min = round(cap_gain_min,2);
cap_gain_max = round(cap_gain_max,2);
% capital loss statistics
cap_loss_mean = round(cap_loss_mean,2);
cap_loss_std = round(cap_loss_std,2);
cap_loss_min = round(cap_loss_min,2);
cap_loss_max = round(cap_loss_max,2);
% age statistics
age_mean = round(age_mean,2);
age_std = round(age_std,2);
age_min = round(age_min,2);
age_max = round(age_max,2);
% hours_per_week_statistics
hpw_mean = round(hpw_mean,2);
hpw_std = round(hpw_std,2);
hpw_min = round(hpw_min,2);
hpw_max = round(hpw_max,2);

cg_stats = table(cap_gain_mean',cap_gain_std','VariableNames',{'mean cap. gains','std cap. gains'},'RowNames',adult_numeric_cat(:,1));
cl_stats = table(cap_loss_mean',cap_loss_std','VariableNames',{'mean cap. loss','std cap. loss'},'RowNames',adult_numeric_cat(:,1));
age_stats = table(age_mean',age_std','VariableNames',{'mean age','std age'},'RowNames',adult_numeric_cat(:,1));
hpw_stats = table(hpw_mean',hpw_std','VariableNames',{'mean h.p.w.','std h.p.w.'},'RowNames',adult_numeric_cat(:,1));

num_stats = [cg_stats, cl_stats, age_stats, hpw_stats];
cap_stats = [cg_stats, cl_stats];
dem_stats = [age_stats, hpw_stats];

num_workclass = num_stats(1:7,:);
num_education = num_stats(8:23,:);
num_mar_stat = num_stats(24:30,:);
num_occupation = num_stats(31:44,:);
num_relationship = num_stats(45:50,:);
num_race = num_stats(51:55,:);
num_gender = num_stats(56:57,:);
num_nc = num_stats(58:59,:);

cap_workclass = cap_stats(1:7,:);
cap_education = cap_stats(8:23,:);
cap_mar_stat = cap_stats(24:30,:);
cap_occupation = cap_stats(31:44,:);
cap_relationship = cap_stats(45:50,:);
cap_race = cap_stats(51:55,:);
cap_gender = cap_stats(56:57,:);
cap_nc = cap_stats(58:59,:);

dem_workclass = dem_stats(1:7,:);
dem_education = dem_stats(8:23,:);
dem_mar_stat = dem_stats(24:30,:);
dem_occupation = dem_stats(31:44,:);
dem_relationship = dem_stats(45:50,:);
dem_race = dem_stats(51:55,:);
dem_gender = dem_stats(56:57,:);
dem_nc = dem_stats(58:59,:);


table2latex(num_workclass,'latextables/num_workclass')
table2latex(num_education,'latextables/num_education')
table2latex(num_mar_stat,'latextables/num_mar_stat')
table2latex(num_occupation,'latextables/num_occupation')
table2latex(num_relationship,'latextables/num_relationship')
table2latex(num_race,'latextables/num_race')
table2latex(num_gender,'latextables/num_gender')
table2latex(num_nc,'latextables/num_nc')

table2latex(cap_workclass,'latextables/cap_workclass')
table2latex(cap_education,'latextables/cap_education')
table2latex(cap_mar_stat,'latextables/cap_mar_stat')
table2latex(cap_occupation,'latextables/cap_occupation')
table2latex(cap_relationship,'latextables/cap_relationship')
table2latex(cap_race,'latextables/cap_race')
table2latex(cap_gender,'latextables/cap_gender')
table2latex(cap_nc,'latextables/cap_nc')

table2latex(dem_workclass,'latextables/dem_workclass')
table2latex(dem_education,'latextables/dem_education')
table2latex(dem_mar_stat,'latextables/dem_mar_stat')
table2latex(dem_occupation,'latextables/dem_occupation')
table2latex(dem_relationship,'latextables/dem_relationship')
table2latex(dem_race,'latextables/dem_race')
table2latex(dem_gender,'latextables/dem_gender')
table2latex(dem_nc,'latextables/dem_nc')

cg_workclass = cg_stats(1:7,:);
cg_education = cg_stats(8:23,:);
cg_mar_stat = cg_stats(24:30,:);
cg_occupation = cg_stats(31:44,:);
cg_relationship = cg_stats(45:50,:);
cg_race = cg_stats(51:55,:);
cg_gender = cg_stats(56:57,:);
cg_nc = cg_stats(58:59,:);

table2latex(cg_workclass,'latextables/cg_workclass')
table2latex(cg_education,'latextables/cg_education')
table2latex(cg_mar_stat,'latextables/cg_mar_stat')
table2latex(cg_occupation,'latextables/cg_occupation')
table2latex(cg_relationship,'latextables/cg_relationship')
table2latex(cg_race,'latextables/cg_race')
table2latex(cg_gender,'latextables/cg_gender')
table2latex(cg_nc,'latextables/cg_nc')

cl_workclass = cl_stats(1:7,:);
cl_education = cl_stats(8:23,:);
cl_mar_stat = cl_stats(24:30,:);
cl_occupation = cl_stats(31:44,:);
cl_relationship = cl_stats(45:50,:);
cl_race = cl_stats(51:55,:);
cl_gender = cl_stats(56:57,:);
cl_nc = cl_stats(58:59,:);


table2latex(cl_workclass,'latextables/cl_workclass')
table2latex(cl_education,'latextables/cl_education')
table2latex(cl_mar_stat,'latextables/cl_mar_stat')
table2latex(cl_occupation,'latextables/cl_occupation')
table2latex(cl_relationship,'latextables/cl_relationship')
table2latex(cl_race,'latextables/cl_race')
table2latex(cl_gender,'latextables/cl_gender')
table2latex(cl_nc,'latextables/cl_nc')

age_workclass = age_stats(1:7,:);
age_education = age_stats(8:23,:);
age_mar_stat = age_stats(24:30,:);
age_occupation = age_stats(31:44,:);
age_relationship = age_stats(45:50,:);
age_race = age_stats(51:55,:);
age_gender = age_stats(56:57,:);
age_nc = age_stats(58:59,:);


table2latex(age_workclass,'latextables/age_workclass')
table2latex(age_education,'latextables/age_education')
table2latex(age_mar_stat,'latextables/age_mar_stat')
table2latex(age_occupation,'latextables/age_occupation')
table2latex(age_relationship,'latextables/age_relationship')
table2latex(age_race,'latextables/age_race')
table2latex(age_gender,'latextables/age_gender')
table2latex(age_nc,'latextables/age_nc')

hpw_workclass = hpw_stats(1:7,:);
hpw_education = hpw_stats(8:23,:);
hpw_mar_stat = hpw_stats(24:30,:);
hpw_occupation = hpw_stats(31:44,:);
hpw_relationship = hpw_stats(45:50,:);
hpw_race = hpw_stats(51:55,:);
hpw_gender = hpw_stats(56:57,:);
hpw_nc = hpw_stats(58:59,:);

table2latex(hpw_workclass,'latextables/hpw_workclass')
table2latex(hpw_education,'latextables/hpw_education')
table2latex(hpw_mar_stat,'latextables/hpw_mar_stat')
table2latex(hpw_occupation,'latextables/hpw_occupation')
table2latex(hpw_relationship,'latextables/hpw_relationship')
table2latex(hpw_race,'latextables/hpw_race')
table2latex(hpw_gender,'latextables/hpw_gender')
table2latex(hpw_nc,'latextables/hpw_nc')
% save every table to latex

% run the acc metrics for each of them
% see the mean and variance of the numeric features for each case
% put them in graphs and get to a conclusion

% PUT IT ALL IN ONE TABLE FOR THE ACCURACY METRICS 
% STUDENT TTEST IS GOOD FOR ALL CASES

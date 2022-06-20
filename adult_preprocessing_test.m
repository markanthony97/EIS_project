%% import table
%clear
adult = readtable('adult.csv');
adult_test = readtable('adult_test.csv');
adult = [adult; adult_test];
adult.Properties.VariableNames = {'age', 'workclass', 'fnlwgt',...
    'education', 'education_num', 'marital_status', 'occupation', ...
    'relationship', 'race', 'sex', 'capital_gain', 'capital_loss', ...
    'hours_per_week', 'native_country', 'result'};
adult = removevars(adult,{'fnlwgt'}); %not important, it indicates how widespread is this kind of person in the census
adult = removevars(adult,{'education_num'}); %not important either, education is already encoded

adult = standardizeMissing(adult, {'?'}); %find and replace ? with missing value

cate_feats = {'workclass', 'education', 'marital_status', 'occupation', ...
    'relationship', 'race', 'sex','native_country','result'};
num_feats = {'age','capital_gain',...
    'capital_loss','hours_per_week'};
adult_rencoded = table();

for j=1:length(num_feats)
    adult_rencoded(:,num_feats{j}) = array2table(...
        table2array(adult(:,num_feats{j})),'VariableNames',{num_feats{j}});
end

for j=1:length(cate_feats)
    adult_rencoded(:,cate_feats{j}) = array2table(categorical(...
        table2array(adult(:,cate_feats{j}))),'VariableNames',{cate_feats{j}});
end

adult_rencoded.result = (grp2idx(adult_rencoded.result) - 1) > 0 ;

findnans=ismissing(adult_rencoded);
find_rows=sum(findnans,2)>0; %find rows with missing values
adult_full = adult_rencoded(~find_rows,:);

%% adult separate numerical and categorical values
cate_feats = cate_feats(1:8);
adult_cat_ohe = table();

for i=1:8
    temp = onehotencode(adult_full(:,cate_feats{i}));
    adult_cat_ohe = [adult_cat_ohe,temp];
end

%% first original trial
% make a difference of capital gains and capital losses
% obtain the PCA
% fcm

adult_numeric = adult_full(:, num_feats);
adult_numeric.net_capital = adult_numeric.capital_gain - adult_numeric.capital_loss;
neo_num_feats = {'age','net_capital','hours_per_week'};
neo_adult_numeric = adult_numeric(:,neo_num_feats);
[~,score,~,tsquared,explained] = pca(table2array(neo_adult_numeric));
fcm_silhouette_vals = zeros(3,11);
for num_c=2:12
    
    [centers,U] = fcm(score,num_c);
    [~,imaxU] = max(U);

    figure(num_c-1)
    subplot(1,2,1)
    hold on
    colors_1 = [0 0 1; 1 0 0; 0 1 0; 1 0.5 0; 0 0.5 1; 0 1 0.5]; 
    colors_2 = [0 0 0.4; 0.4 0 0; 0 0.4 0; 0.4 0.2 0; 0 0.2 0.4; 0 0.4 0.2];
    clusterresperc = zeros(1,num_c);
    clust_names = cell(1,num_c);
    
    for i=1:num_c
         colortrial = rand(1,1,3);
         scatter3(score(imaxU==i,1),score(imaxU==i,2),score(imaxU==i,3),'MarkerEdgeColor',colortrial)
         clusterresperc(i) = sum(adult_full.result(imaxU==i))/length(adult_full.result(imaxU==i))*100;
         display(strcat("La proporzione di campioni appartenenti al cluster ",...
                 num2str(i),  " che hanno valore di target pari a >=50K è ",...
                 num2str(clusterresperc(i))," %"))
         clust_names(i) = {strcat('cluster ',num2str(i))};
    end

    for i=1:num_c
        scatter3(centers(i,1),centers(i,2),centers(i,3),100,'+','LineWidth',1,'MarkerFaceColor',colortrial,'MarkerEdgeColor','k')
    end
    title(strcat('Clustering of ',num2str(num_c),' clusters'))
    grid
    legend(clust_names)
    hold off 
    subplot(1,2,2)
    title('silhouette values')
    silhouette(score,imaxU)
    fcm_silhouette_vals(1,num_c)=mean(silhouette(score,imaxU));
end
%% with the steps
% 1 normalization
% 2 pca
% 3 fcm
neo_adult_numeric_2 = array2table(zscore(table2array(neo_adult_numeric)),'VariableNames', neo_num_feats);
[~,score_2,~,tsquared_2,explained_2] = pca(table2array(neo_adult_numeric_2));
%%
for num_c=2:12
    
    [centers_2,U_2] = fcm(score_2,num_c);
    [~,imaxU_2] = max(U_2);

    figure(num_c-1)
    subplot(1,2,1)
    hold on
    clust_names = cell(1,num_c);
    
    for i=1:num_c
         colortrial = rand(1,1,3);
         scatter3(score_2(imaxU_2==i,1),score_2(imaxU_2==i,2),score_2(imaxU_2==i,3),'MarkerEdgeColor',colortrial)
         clusterresperc = sum(adult_full.result(imaxU_2==i))/length(adult_full.result(imaxU_2==i))*100;
         display(strcat("La proporzione di campioni appartenenti al cluster ",...
                 num2str(i),  " che hanno valore di target pari a >=50K è ",...
                 num2str(clusterresperc)," %"))
         clust_names(i) = {strcat('cluster ',num2str(i))};
    end

    for i=1:num_c
        scatter3(centers_2(i,1),centers_2(i,2),centers_2(i,3),100,'+','LineWidth',1,'MarkerFaceColor',colortrial,'MarkerEdgeColor','k')
    end
    title(strcat('Clustering of ',num2str(num_c),' clusters'))
    grid
    legend(clust_names)
    hold off 
    subplot(1,2,2)
    title('silhouette values')
    silhouette(score_2,imaxU_2)
    
    fcm_silhouette_vals(2,num_c)=mean(silhouette(score_2,imaxU_2));
end
   

%% 
adult_numeric_3 = zscore(table2array(neo_adult_numeric));
for num_c=2:12
    
    [centers_2,U_2] = fcm(adult_numeric_3,num_c);
    [~,imaxU_3] = max(U_2);

    figure(num_c-1)
    subplot(1,2,1)
    hold on
    clust_names = cell(1,num_c);
    
    for i=1:num_c
         colortrial = rand(1,1,3);
         scatter3(adult_numeric_3(imaxU_3==i,1),adult_numeric_3(imaxU_3==i,2),adult_numeric_3(imaxU_3==i,3),'MarkerEdgeColor',colortrial)
         clusterresperc = sum(adult_full.result(imaxU_3==i))/length(adult_full.result(imaxU_3==i))*100;
         display(strcat("La proporzione di campioni appartenenti al cluster ",...
                 num2str(i),  " che hanno valore di target pari a >=50K è ",...
                 num2str(clusterresperc)," %"))
         clust_names(i) = {strcat('cluster ',num2str(i))};
    end

    for i=1:num_c
        scatter3(centers_2(i,1),centers_2(i,2),centers_2(i,3),100,'+','LineWidth',1,'MarkerFaceColor',colortrial,'MarkerEdgeColor','k')
    end
    title(strcat('Clustering of ',num2str(num_c),' clusters'))
    grid
    legend(clust_names)
    hold off 
    subplot(1,2,2)
    title('silhouette values')
    silhouette(adult_numeric_3,imaxU_3)
    
    fcm_silhouette_vals(3,num_c)=mean(silhouette(adult_numeric_3,imaxU_3));
end

%% CONFERMATO CHE LA NORMALIZZAZIONE PEGGIORA SOLO LE COSE
 
figure(4)
hold on
cluster = zeros(1,num_c);
clust_names = cell(1,num_c);
for i=1:num_c
     scatter3(score_2(imaxU_2==i,1),score_2(imaxU_2==i,2),score_2(imaxU_2==i,3),'MarkeredgeColor',colors_1(i,:))
     %scatter3(centers(i,1),centers(i,2),centers(i,3),100,'+','LineWidth',3,'MarkeredgeColor',colors_2(i,:))
     cluster = sum(adult_full.result(imaxU_2==i))/length(adult_full.result(imaxU_2==i))*100;
     display(strcat("La proporzione di campioni appartenenti al cluster ",...
             num2str(i),  " che hanno valore di target pari a >=50K è ",...
             num2str(cluster)," %"))
     clust_names(i) = {strcat('cluster ',num2str(i))};
end
grid
legend(clust_names)
for i=1:num_c
    scatter3(centers(i,1),centers(i,2),centers(i,3),100,'+','LineWidth',3,'MarkeredgeColor',colors_2(i,:))
end
hold off 

% scatter3(score(imaxU==1,1),score(imaxU==1,2),score(imaxU==1,3),'ob')
% hold on
% scatter3(score(imaxU==2,1),score(imaxU==2,2),score(imaxU==2,3),'or')
% scatter3(score(imaxU==3,1),score(imaxU==3,2),score(imaxU==3,3),'og')
% scatter3(score(imaxU==4,1),score(imaxU==4,2),score(imaxU==4,3),'o','MarkeredgeColor',[1 0.5 0])
% scatter3(centers(1,1),centers(1,2),centers(1,3),100,'+','LineWidth',3,'MarkeredgeColor',[0 0 0.4])
% scatter3(centers(2,1),centers(2,2),centers(2,3),100,'+','LineWidth',3,'MarkeredgeColor',[0.4 0 0])
% scatter3(centers(3,1),centers(3,2),centers(3,3),100,'+','LineWidth',3,'MarkeredgeColor',[0 0.4 0])
% scatter3(centers(4,1),centers(4,2),centers(4,3),100,'+','LineWidth',3,'MarkeredgeColor',[0.8 0.5 0])
% legend('cluster 1','cluster 2','cluster 3','cluster 4')
% hold off

%% 
% use the result variable to consider the clusters valid or the whole
% dataset?
cluster(1) = sum(adult_full.result(imaxU==1))/length(adult_full.result(imaxU==1))*100;
cluster(2) = sum(adult_full.result(imaxU==2))/length(adult_full.result(imaxU==2))*100;
% cluster(3) = sum(adult_full.result(imaxU==3))/length(adult_full.result(imaxU==3))*100;
% cluster(4) = sum(adult_full.result(imaxU==4))/length(adult_full.result(imaxU==4))*100;
cluster = round(cluster,2);


% clusters = categorical(imaxU',[1,2,3,4],{'c_1','c_2','c_3','c_4'});
clusters = categorical(imaxU',[1,2],{'c_1','c_2'});
% cluster_ohe = onehotencode(clusters,2);


% adult_full_ohe = [cluster_ohe, table2array(adult_cat_ohe)];


%% Find exceptions and maximum to begin clustering
eq_matrix = adult_full_ohe*adult_full_ohe';

lone_wolves = sum(eq_matrix==0);
[whatisthis,black_sheep] = maxk(lone_wolves,10);

adult_full(black_sheep,:)
whatisthis

%%

eq_mat = eq_matrix.*~eye(size(eq_matrix));
max_found = max(eq_mat,[],'all'); % this shifts constantly

max_eq = eq_matrix==max_found;


%questo è il nuovo input alla iterazione successiva, ripeti pca e fcm su di
%questo per escludere i cluster

adult_new=adult_full(~max_eq(:,I),:); 

% TODO generalizza la procedura per poter escludere via via le cifre e
% rifare il calcolo, alla fine ci si dovrebbe ritrovare con un certo numero
% di raggruppamenti


%% wrapper functions


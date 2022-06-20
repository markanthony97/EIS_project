cate_feats = {'workclass', 'education', 'marital_status', 'occupation', ...
    'relationship', 'race', 'sex', 'native_country'};
num_feats = {'capital_gain','capital_loss','age',...
    'hours_per_week'};
load grouped_clusters
unique_countries = unique(grouped_clusters.native_country);
n_of_occurrences = zeros(1,length(unique_countries));
for i=1:length(unique_countries)
    n_of_occurrences(i) = sum(grouped_clusters.native_country==unique_countries(i));
end
countries = table(unique_countries, n_of_occurrences','VariableNames',{'countries','n_of_records'});
countries = sortrows(countries,'n_of_records','ascend');
grouped_clusters.native_country(grouped_clusters.native_country~='United-States') = 'Foreign-Country';
lone_clusters.native_country(lone_clusters.native_country~='United-States') = 'Foreign-Country';


[adult_for_nn, result, adult_numeric, adult_categoric] = prepare_adult(grouped_clusters, num_feats, cate_feats);
[adult_test, result_test, lone_numeric, lone_categoric] = prepare_adult(lone_clusters, num_feats, cate_feats);
adult_for_nn_reduced = adult_for_nn(:,sum(adult_for_nn>0,1)>0);
adult_test = adult_test(:,sum(adult_test>0,1)>0);
function [adult_for_nn ,result, adult_numeric, adult_categoric] = prepare_adult(grouped_clusters, num_feats, cate_feats)
    adult_numeric = grouped_clusters(:,num_feats);
    adult_categoric = grouped_clusters(:,cate_feats);
    result = grouped_clusters.result;

    adult_categoric_ohe = table();
    for i=1:8
        temp = onehotencode(adult_categoric(:,cate_feats{i}));
        adult_categoric_ohe = [adult_categoric_ohe, temp];
    end

    % adult_numeric_normalized = array2table(zscore(table2array(adult_numeric)),'VariableNames', num_feats);
    adult_for_nn = [adult_numeric,adult_categoric_ohe];
    adult_for_nn=table2array(adult_for_nn);
end
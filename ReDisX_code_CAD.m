%%% clustering to the cobat output
clear all;close all;
addpath '/home/yip/Documents/TVRFtest (1)/TVRFtest'
addpath '/home/yip/Documents/Reclassification_customcdf/SupervisedPCA-master'
rng('default')



CAD_data = readtable('Top5000CAD_data.csv');
CAD_expression=table2array(CAD_data(:,1:end));
CAD_sample_id=CAD_data.Properties.VariableNames(1:end);
CAD_gene_name=readtable('Top5000CAD_genename.csv');
CAD_gene_name.Properties.VariableNames={'SYMBOL'};



CAD_phenodata = readtable('phenotype_CAD_data.csv');

CAD_disease_type=categorical(CAD_phenodata.disease_type);
CAD_control_=find(CAD_disease_type == "Control");





%%% only control
control_expression=CAD_expression(:,CAD_control_);
control_sample_id=CAD_sample_id(:,CAD_control_);
control_phenodata=CAD_phenodata(CAD_control_,:);
control_disease_type=CAD_disease_type(CAD_control_);


CAD_disease_type(CAD_control_)=[];
CAD_expression(:,CAD_control_)=[];
CAD_sample_id(:,CAD_control_)=[];
CAD_phenodata(CAD_control_,:)=[];






%%% all_in_one

all_expression=[CAD_expression control_expression];
all_sample_id=[CAD_sample_id control_sample_id];
all_phenodata=[CAD_phenodata;control_phenodata ];
all_disease_type=[CAD_disease_type;control_disease_type];


Normalized_CAD_exp = normalize(CAD_expression,2);
T_eva=[];
predicted_by_CMF_eva=[];
C=[];
max_clus=20;
%%% evaluate number of cluster by ReDisX
for i=1:max_clus
    T_eva{i} = clusterdata(transpose(Normalized_CAD_exp),'distance','correlation','linkage','ward','Maxclust',i);
    [GC,GR] = groupcounts(T_eva{i});
    CMF_numberoflabel(i)=min(GC);
    predicted_by_CMF_eva{i}=maxflow_main(transpose(Normalized_CAD_exp),i,T_eva{i},1,CMF_numberoflabel(i),'log_region_force_2');
    for j=1:length(CAD_sample_id)
     [m,final_clus_tem(i,j)]=max(predicted_by_CMF_eva{i}(j,:));
    end
    C(i) = length(unique(final_clus_tem(i,:)));

end
[GC,GR] = groupcounts(transpose(C));

[maxnum,maxin]=max(GC(2:end));
true_optimal_cluster=GR(maxin+1);
optimal_cluster=find(C==5,1);


figure_eva=figure;
bar(GR(2:end),GC(2:end)/max_clus,'DisplayName','GC')
ylabel('frequenct')
xlabel('number of cluster')
title('cluster evaluation by ReDisX')
saveas(figure_eva,'result_v2/ReDisX_GSE59867_eva_clus.png')





T = T_eva{optimal_cluster};
predicted_by_CMF=predicted_by_CMF_eva{optimal_cluster};
final_clus=zeros(length(CAD_sample_id),1);
for i=1:length(CAD_sample_id)
[m,final_clus(i)]=max(predicted_by_CMF(i,:));
end
index_six = find(final_clus==6);
final_clus(final_clus==6)=5;


Y = tsne(transpose(Normalized_CAD_exp),'NumPCAComponents',50);

figure1=figure;
gscatter(Y(:,1),Y(:,2),final_clus);
title('2-D Embedding ReDisX result')
saveas(figure1,'result_v2/ReDisX_GSE59867_result.png')







final_table_result=[table(transpose(CAD_sample_id),final_clus, CAD_disease_type) CAD_phenodata(:,3)];

writetable(final_table_result,'matlab_GSE59867_ReDisX/CAD_clustering_result.csv')
writematrix(CAD_expression,'matlab_GSE59867_ReDisX/CAD_expression.csv')




final_table_withcontrol_result=[table(transpose(all_sample_id),[final_clus;(true_optimal_cluster+1)*ones(length(CAD_control_),1)] , all_disease_type) all_phenodata(:,3)];
writetable(final_table_withcontrol_result,'matlab_GSE59867_ReDisX/CAD_clustering_result_pluscontrol.csv')
writematrix(all_expression,'matlab_GSE59867_ReDisX/CAD_expression_final_pluscontrol.csv')

save matlab_GSE59867_ReDisX/clustering_result_CAD.mat


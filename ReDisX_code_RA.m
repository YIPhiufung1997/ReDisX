%%% clustering to the R output
clear all;close all;
addpath '/home/yip/Documents/TVRFtest (1)/TVRFtest'

rng('default')



RA_data = readtable('Top5000RA_data.csv');
RA_expression=table2array(RA_data(:,1:end));
RA_sample_id=RA_data.Properties.VariableNames(1:end);
RA_gene_name=readtable('Top5000RA_genename.csv');
RA_gene_name.Properties.VariableNames={'SYMBOL'};



RA_phenodata = readtable('phenotype_RA_data.csv');

RA_disease_type=categorical(RA_phenodata.disease_type);
RA_control_=find(RA_disease_type == "Control");





%%% only control
control_expression=RA_expression(:,RA_control_);
control_sample_id=RA_sample_id(:,RA_control_);
control_phenodata=RA_phenodata(RA_control_,:);
control_disease_type=RA_disease_type(RA_control_);


RA_disease_type(RA_control_)=[];
RA_expression(:,RA_control_)=[];
RA_sample_id(:,RA_control_)=[];
RA_phenodata(RA_control_,:)=[];






%%% all_in_one

all_expression=[RA_expression control_expression];
all_sample_id=[RA_sample_id control_sample_id];
all_phenodata=[RA_phenodata;control_phenodata ];
all_disease_type=[RA_disease_type;control_disease_type];


Normalized_RA_exp = normalize(RA_expression,2);


%%% evaluate number of cluster by CMF
T_eva=[];
predicted_by_CMF_eva=[];
C=[];
max_clus=20;
CMF_numberoflabel=[];
for i=1:max_clus
    T_eva{i} = clusterdata(transpose(Normalized_RA_exp),'distance','correlation','linkage','ward','Maxclust',i);
    [GC,GR] = groupcounts(T_eva{i});
    CMF_numberoflabel(i)=min(GC);
    predicted_by_CMF_eva{i}=maxflow_main(transpose(Normalized_RA_exp),i,T_eva{i},1,CMF_numberoflabel(i),'log_region_force_2');
    for j=1:length(RA_sample_id)
     [m,final_clus_tem(i,j)]=max(predicted_by_CMF_eva{i}(j,:));
    end
    C(i) = length(unique(final_clus_tem(i,:)));

end
[GC,GR] = groupcounts(transpose(C));

[maxnum,maxin]=max(GC(2:end));
optimal_cluster=GR(maxin+1);

figure_eva=figure;
bar(GR(2:end),GC(2:end)/max_clus,'DisplayName','GC')
ylabel('frequenct')
xlabel('number of cluster')
title('cluster evaluation by ReDisX')
saveas(figure_eva,'result_v2/ReDisX_GSE93272_eva_clus.png')






T = T_eva{optimal_cluster};
predicted_by_CMF=predicted_by_CMF_eva{optimal_cluster};
final_clus=zeros(length(RA_sample_id),1);
for kk=1:length(RA_sample_id)
[m,final_clus(kk)]=max(predicted_by_CMF(kk,:));
end



Y = tsne(transpose(Normalized_RA_exp),'NumPCAComponents',50);

figure1=figure;
gscatter(Y(:,1),Y(:,2),final_clus);
title('2-D Embedding ReDisX result')
saveas(figure1,'result_v2/ReDisX_GSE93272_result.png')







final_table_result=[table(transpose(RA_sample_id),final_clus, RA_disease_type) RA_phenodata(:,3)];

writetable(final_table_result,'matlab_GSE93272_ReDisX/RA_clustering_result.csv')
writematrix(RA_expression,'matlab_GSE93272_ReDisX/RA_expression.csv')




final_table_withcontrol_result=[table(transpose(all_sample_id),[final_clus;(optimal_cluster+1)*ones(length(RA_control_),1)] , all_disease_type) all_phenodata(:,3)];
writetable(final_table_withcontrol_result,'matlab_GSE93272_ReDisX/RA_clustering_result_pluscontrol.csv')
writematrix(all_expression,'matlab_GSE93272_ReDisX/RA_expression_final_pluscontrol.csv')

save matlab_GSE93272_ReDisX/clustering_result_RA.mat


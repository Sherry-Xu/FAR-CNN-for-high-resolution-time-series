node_list  = readtable('../data_new/nodes_list.csv','ReadVariableNames',false,'ReadRowNames',false);
node_list2 = table2cell(node_list);
for i  = 1:length(node_list2)
    index = node_list2(i);
    index = char(index);
    path = strcat('../data_reshape/',index,'.csv');
    data = csvread(path);
    FAR(data,strcat('data_FAR/',index,'_FARout.mat'),strcat('data_FAR/',index,'_res.csv'),strcat('data_FAR/',index,'_yhat.csv'),11,28)
end

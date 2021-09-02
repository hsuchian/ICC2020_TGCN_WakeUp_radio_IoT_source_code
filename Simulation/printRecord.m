% Driverz
average_num = 10;

for fi = 1:0.5:1
    fileID = fopen(['exp', num2str(fi),'.m'],'a');
    for i = 121:200
        [SuccessProb, record, estimated_success_prob] = SimHash1(i, round(i*fi), -1000, average_num);
        fprintf(fileID, '%d ', estimated_success_prob);
    end
    fclose(fileID);
end

%fileID = fopen('exp.m','w');
%fprintf(fileID, '%d ', estimated_success_prob);
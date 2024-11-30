
function [r2,rmse] = cal_r(X,Y)

if size(X)~=size(Y)
    error 'X size is unequal to Y';
end

r2_n = [];rmse_n = [];
for n = 1:size(X,1)
    [r2_n(n),rmse_n(n)] = rsquare(X(n,:),Y(n,:));
end

% r2=[mean(r2_n);std(r2_n)];
% rmse=[mean(rmse_n);std(rmse_n)];

r2=mean(r2_n);
rmse=mean(rmse_n);

end
% 权重剪枝 以权重数值为阈值进行剪枝
function w_pruning=weight_pruning_threshold(w,beta)
% s为一个百分数
W=sort(w,'descend');
loc=find(abs(w)>0);
w0=length(loc);
% beta为阈值
% 进行剪枝操作
for i=1:length(w)
    if w(i)<beta
        w(i)=0;
    end
end
w_pruning=w;

% 重新训练
loc=find(w==0);
% w = zeros(tapslen_1+tapslen_2*(tapslen_2+1)/2 + tapslen_3*(tapslen_3+1)*(tapslen_3+2)/6 +...
%     fblen_1+fblen_2*(fblen_2+1)/2 + fblen_3*(fblen_3+1)*(fblen_3+2)/6 ,1);
index_1=[];
index_2=[];
index_3=[];
for i=1:length(loc)
    L=loc(i);
    % 待丢弃的索引
    if L<tapslen_1
        a=L;
        index_1=[index_1 a];
    elseif L<tapslen_1+tapslen_2*(tapslen_2+1)/2
        b=L-tapslen_1;
        index_2=[index_2 b];
    else
        c=L-(tapslen_1+tapslen_2*(tapslen_2+1)/2);
        index_3=[index_3 c];
    end

end
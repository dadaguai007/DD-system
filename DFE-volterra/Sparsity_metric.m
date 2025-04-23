% Volterra Sparsity metric
% 计算Volterra的稀疏度
function phi=Sparsity_metric(w,taps_list)
NC=[];
for i=1:length(taps_list)

    c1=factorial(taps_list(i)+i-1);
    c2=factorial(taps_list(i)-1);
    c3=factorial(i);
    Mu=c1./(c2.*c3);
    NC=[NC Mu];

end
NC=sum(NC);
loc=abs(w)>0;
w1=w(loc);
phi=1-(w1)./NC;

end

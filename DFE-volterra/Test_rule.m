% Test rule
clc;clear;
x=0:0.01:1;
A=[5,5,10];B=[0.5,0.2,0.5];
% A=5;B=0.2;
for i=1:length(A)
    a=A(i);
    b=B(i);
x1=1-exp(-a*(x./b-1));
x2=1+exp(-a*(x./b-1));
f(i,:)=1/2*(x1./x2+1);
end

f_rule1=ones(1,length(x));
f_rule1(1:floor(length(x)/2)+1)=0;

f_rule2=x;
figure;hold on
plot(x,f(1,:),'LineWidth',2)
plot(x,f(2,:),'LineWidth',2)
plot(x,f(3,:),'LineWidth',2)
plot(x,f_rule1,'LineWidth',2)
plot(x,f_rule2,'LineWidth',2)
legend('1','2','3','4','5')
Plotter('function comparison for different rules.','x','f(x)');
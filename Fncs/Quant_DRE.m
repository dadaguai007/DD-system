function [sigout] = Quant_DRE(sigin,h,Qbit_num)
%%   根据论文《Low-Resolution Digital Pre-Compensation Enabled by Digital Resolution Enhancer》逻辑编写
%%   author：cc
%%   PNOB默认大于1
if size(sigin,1) > 1
    sigin = sigin.';
end

sigini=real(sigin);
% Qbit_num=floor(2^pnob);
% Qbit_num=floor(2^pnob);
I_max = max(abs(sigini));
I_min = -I_max;




step = (I_max-I_min)/Qbit_num;
for ii = 1:Qbit_num-1
    ipartition(1,ii) = I_min+step*ii;
end
for jj = 1:Qbit_num
    icodebook(1,jj)  = I_min + step/2 + step*(jj-1);
end
[Iorder,iroundsig]=quantiz(sigini,ipartition,icodebook);

Iorder = Iorder+1;

I_error = iroundsig-sigini;

%% ...

n=length(icodebook);
m=length(h);
nvtb=length(sigin);
I_error = [I_error,I_error(1:m)];
Iorder = [Iorder,Iorder(1:m)];
Qro1 = zeros(3,m);
Qsq1 = zeros(3,m);
Qro2 = zeros(2,m);
Qsq2 = zeros(2,m);
for ii = m+1:nvtb+m;
    if Iorder(ii) ~= 1 && Iorder(ii) ~= Qbit_num  %% 判断是否是边帧，因为边帧的话只有两种可能性
        Qro1(1,:) = I_error(ii-m+1:ii)-[zeros(1,m-1),step];
        Qro1(2,:) = I_error(ii-m+1:ii);
        Qro1(3,:) = I_error(ii-m+1:ii)+[zeros(1,m-1),step];
        Qsq1(1,:) = conv(Qro1(1,:),h,'same');
        Qsq1(2,:) = conv(Qro1(2,:),h,'same');
        Qsq1(3,:) = conv(Qro1(3,:),h,'same');
        [~,num] = min(sum(abs(Qsq1),2));
        I_error(ii) = Qro1(num,m);
    else
        Qro2(1,:) = I_error(ii-m+1:ii)-[zeros(1,m-1),step].*sign(Iorder(ii)-n/2); %% 减少if等判断语句，加快程序运行速度
        % 如果是最低位，则考虑最低和次低位，如果是最高位，则考虑最高和次高位
        Qro2(2,:) = I_error(ii-m+1:ii);
        Qsq2(1,:) = conv(Qro2(1,:),h,'same');
        Qsq2(2,:) = conv(Qro2(2,:),h,'same');
        [~,num] = min(sum(abs(Qsq2),2));
        I_error(ii) = Qro2(num,m);
    end
end

I_error1 = [I_error(end-m+1:end),I_error(m+1:nvtb)];

sigini1 = I_error1+sigini;

sigout = sigini1;
% scatterplot(sigout)


if size(sigout,2) > 1
    sigout = sigout.';
end

end
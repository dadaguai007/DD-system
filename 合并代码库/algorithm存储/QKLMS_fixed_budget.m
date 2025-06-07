function[y,e,time] = QKLMS_fixed_budget(nsig,d,step,tap_num,delay,sps,db)
if nargin < 6
    sps = 2;
end
a = 0.2; %tap_num越多，a越小
epsilon = 1;
% db = 2000;
nsig = nsig(:)';
d = d(:)';
u_all = buffer(nsig,tap_num,tap_num-sps,"nodelay");
%%迭代滤波%%
e=zeros(1,size(u_all,2));%误差矩阵
y=zeros(1,size(u_all,2));%滤波输出
z(1) = step*d(1);
c(:,1) = u_all(:,1);
tic
m = 2;
for i=2:size(u_all,2)%开始迭代更新参数
%     while m <= db
        k_gauss = exp(-a*sum((c-repmat(u_all(:,i),1,m-1)).^2));
        y(i-1)=z*k_gauss';%输出
        e(i-1)=d(i+delay)-y(i-1);%误差
        [dis,idx] = min(sqrt(sum((c-repmat(u_all(:,i),1,m-1)).^2)));
        if dis <= epsilon
            z(idx) = z(idx)+step*e(i-1);
        else
            c = [c u_all(:,i)];
            z(m) = step*e(i-1);
            m = m+1;
        end
        if m > db
            break;
        end
%     end
end
for i = i+1:size(u_all,2)
    k_gauss = exp(-a*sum((c-repmat(u_all(:,i),1,m-1)).^2));
    y(i-1)=z*k_gauss';%输出
    e(i-1)=d(i+delay)-y(i-1);%误差
    [dis,idx] = min(sqrt(sum((c-repmat(u_all(:,i),1,m-1)).^2)));
    if dis <= epsilon
        z(idx) = z(idx)+step*e(i-1);
    else
        c(:,idx) = u_all(:,i);
        z(idx) = z(idx)+step*e(i-1);

% c = [c(:,(2:end)) u_all(:,i)];
% z = [z(2:end) step*e(i-1)];
    end
end
time = toc;
if size(c,2)<db
    warning('Space budget not fully utilized');
end
% plot(n);
end
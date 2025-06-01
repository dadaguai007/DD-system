function plotEquParam(xn,sigRx_E,label,w,en)

% 均衡后信号时域分布
figure;
hold on;
plot(xn(1:1e5),'.')
plot(sigRx_E(1:1e5),'k.')
plot(label(1:1e5),'.')
legend('接收信号','均衡后信号','发送信号')
const=unique(label);
ylim([min(const)-2 max(const)+2])
% 抽头数
figure;
stem(w(:))

% 误差是否收敛
figure;
semilogy(abs(en(1:1e4)).^2)
xlabel("迭代次数")
ylabel("误差")

end
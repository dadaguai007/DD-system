% test comn.MemorylessNonlinear
% 先过放大器的非线性，再过IQ不平衡的非线性
%指定了它的三阶截点（IIP3）为38 dBm，并且假设没有幅相转换（AM-PM Conversion）。
% 增加放大器增益压缩，而不会造成AM-PM转换障碍,后面加上，会有相位偏转
%非线性放大器增益压缩
mnlamp= comm.MemorylessNonlinearity(IIP3=38, AMPMConversion=0);
plot(mnlamp);
%%
M = 16;           % Modulation order
sps = 4;          % Samples per symbol
pindBm = [12 25]; % Input power
gain = 10;        % Amplifier gain
% ReferenceImpedance代表的是输出阻抗。。。
amplifier = comm.MemorylessNonlinearity("Method","Cubic polynomial", ...
    "LinearGain",gain,"AMPMConversion",10,"ReferenceImpedance",50);
pin = 10.^((pindBm-30)/10); % Convert dBm to linear Watts
data = randi([0 M-1],1000,1);
modOut = qammod(data,M,"UnitAveragePower",true)*sqrt(pin*amplifier.ReferenceImpedance);
ampOut = amplifier(modOut);
plot(amplifier);
% 输出信号两个放大功率不一致
figure;hold on;
plot(real(modOut(:,1)));
plot(real(ampOut(:,1)));
%% LUT的放大器模型

pindBm = [-8,0];                % Input power
pin = 10.^((pindBm-30)/10); % power in Watts
paChar = pa_performance_characteristics();
amplifier = comm.MemorylessNonlinearity('Method','Lookup table','Table',paChar,'ReferenceImpedance',50);
modSig = qammod(data,M,'UnitAveragePower',true)*sqrt(pin * amplifier.ReferenceImpedance);
ampSig = amplifier(modSig);
plot(amplifier)
%% 'Saleh model'
pindBm = [-8,0];                % Input power
pin = 10.^((pindBm-30)/10); % power in Watts
saleh = comm.MemorylessNonlinearity('Method','Saleh model','ReferenceImpedance', 50);
modSig1 = qammod(data,M,'UnitAveragePower',true)*sqrt(pin * amplifier.ReferenceImpedance);
ampSig1 = saleh(modSig1);
plot(saleh)



%% Hyper 跟三次多项式通用
pindBm = [-8,0];                % Input power
pin = 10.^((pindBm-30)/10); % power in Watts
hyper = comm.MemorylessNonlinearity('Method','Hyperbolic tangent',...
    IIP3=20, AMPMConversion=0,ReferenceImpedance=50);
modSig2 = qammod(data,M,'UnitAveragePower',true)*sqrt(pin * amplifier.ReferenceImpedance);
ampSig2 = hyper(modSig2);
plot(hyper)

%%
% Ghorbani model
pindBm = [-8,0];                % Input power
pin = 10.^((pindBm-30)/10); % power in Watts
Ghorbani = comm.MemorylessNonlinearity('Method','Ghorbani model','ReferenceImpedance', 50);
modSig3 = qammod(data,M,'UnitAveragePower',true)*sqrt(pin * amplifier.ReferenceImpedance);
ampSig3 = Ghorbani(modSig3);
plot(Ghorbani)
%%
function paChar = pa_performance_characteristics()
HAV08_Table =...
    [-35,60.53,0.01;
    -34,60.53,0.01;
    -33,60.53,0.08;
    -32,60.54,0.08;
    -31,60.55,0.1;
    -30,60.56,0.08;
    -29,60.57,0.14;
    -28,60.59,0.19;
    -27,60.6,0.23;
    -26,60.64,0.21;
    -25,60.69,0.28;
    -24,60.76,0.21;
    -23,60.85,0.12;
    -22,60.97,0.08;
    -21,61.12,-0.13;
    -20,61.31,-0.44;
    -19,61.52,-0.94;
    -18,61.76,-1.59;
    -17,62.01,-2.73;
    -16,62.25,-4.31;
    -15,62.47,-6.85;
    -14,62.56,-9.82;
    -13,62.47,-12.29;
    -12,62.31,-13.82;
    -11,62.2,-15.03;
    -10,62.15,-16.27;
    -9,62,-18.05;
    -8,61.53,-20.21;
    -7,60.93,-23.38;
    -6,60.2,-26.64;
    -5,59.38,-28.75];
paChar = HAV08_Table;
paChar(:,2) = paChar(:,1) + paChar(:,2);
end



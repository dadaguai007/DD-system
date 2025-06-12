classdef DDTx < handle

    % 定义类的属性
    properties
        % 为了方便后续添加变量，采用结构体方式进行变量命名
        TxPHY % 发射机参数
        Nr %无量纲参数
        Button % 开关讯号
        Implementation
    end

    % 定义类的方法
    methods
        % 构造函数
        function obj = DDTx(varargin)
            %初始化类的属性
            if numel(varargin) == 13
                obj.TxPHY.fs                           = varargin{1} ;% 发射信号的采样率
                obj.TxPHY.fb                           = varargin{2} ;% 发射信号的波特率
                obj.TxPHY.order                        = varargin{3} ;% 随机信号的阶数
                obj.TxPHY.prbsOrder                    = varargin{4} ;% prbs码的阶数
                obj.TxPHY.M                            = varargin{5} ;% 调制格式
                obj.TxPHY.sps                          = varargin{6} ;% 每符号采样点
                obj.TxPHY.NSym                         = varargin{7} ;% 码元数目
                obj.Nr.psfShape                        = varargin{8} ;% 脉冲形式
                obj.Nr.psfRollOff                      = varargin{9} ;% 滚降系数
                obj.Nr.psfLength                       = varargin{10};% 影响长度
                obj.Nr.userShape                       = varargin{11};% 用户的成型滤波器
                obj.Button.DataType                    = varargin{12};% 选择模式
                obj.Button.shapingFilter               = varargin{13};% 成型滤波器的生成方式
            end
            obj.Button.delaySignal = 'symbol';          % 采用采样点延迟
            obj.Nr.ncut_index = 0.2*obj.TxPHY.NSym; %解码起始位置
        end

        % 参考信号(解码使用)
        function refOut=createReferenceSignal(obj,ref)
            % 得到用于解码的参考信号
            obj.Implementation.ref = repmat(ref,1,100);
            refOut=obj.Implementation.ref;
        end

        % 创建参考星座图
        function [const,Ksym]=creatReferenceConstellation(obj)
            % Constellation
            z = (0:obj.TxPHY.M-1)';
            % 星座图生成
            const=pammod(z,obj.TxPHY.M);
            Ex       = 1;          % 符号平均能量
            Ksym = modnorm(const, 'avpow', Ex);   % 能量归一化因子
            obj.Implementation.constellation  = Ksym * const; % 调整后的参考星座
        end

        function [decodedData,ber] = PAM_ExecuteDecoding(obj, eq_signal,const)
            % 针对 PAM4 解码，计算误码
            % 参考信号
            ref_seq=obj.Implementation.ref;
            label=decision(ref_seq,const);
            % 参考序列
            label_bit=obj.pam4demod(label);
            % 接收序列
            I=decision(eq_signal,const);
            decodedData=obj.pam4demod(I);
            % 解码
            [ber,num,~] = CalcBER(decodedData(obj.Nr.ncut_index:end),label_bit(obj.Nr.ncut_index:end)); %计算误码率
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        function  [decodedData,ber] = PAM_QuantifyDecoding(obj, eq_signal)
            % 量化区间
            A=[-2 0 2];
            ref_seq=obj.Implementation.ref;
            % 参考序列
            [~,label] = quantiz(ref_seq,A,[-3,-1,1,3]);
            label_bit=obj.pam4demod(label);
            % 接收序列
            if strcmp(obj.Button.Train,'on')
                [~,I] = quantiz(eq_signal,A,[-3,-1,1,3]);
            else
                [~,I] = quantiz(eq_signal,pnorm(A),[-3,-1,1,3]);
            end
            decodedData=obj.pam4demod(I);
            % 解码
            [ber,num,~] = CalcBER(decodedData(obj.Nr.ncut_index:end),label_bit(obj.Nr.ncut_index:end)); %计算误码率
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        % 输出信号
        function [filteredSignal,pam_signal,pulse]=dataOutput(obj)
            % 成型滤波器
            if strcmp(obj.Button.shapingFilter,'system')
                pulse = obj.systemHsqrt();
            elseif strcmp(obj.Button.shapingFilter,'user')
                pulse = obj.getPulse(obj.Nr.userShape);
            end

            %生成bit数
            switch lower(obj.Button.DataType)
                case 'prbs'
                    [~,symbols]=obj.prbs_bits();
                case 'rand'
                    symbols=obj.rand_bits();
            end
            %调制 , 装载信号
            pam_signal=obj.pam(symbols);
            pam_signal=real(pam_signal);
            % 归一化
            pamSignalNorm = pnorm(pam_signal);
            % 脉冲成型 装载成型后的信号
            filteredSignal=obj.applyShapingFilter(pamSignalNorm,pulse);
        end


        function duobSignal=Duobdataout(obj)
            %生成bit数
            symbols=obj.modcode();

            if obj.TxPHY.M > 2
                % 应用格雷编码 ,2,3调换
                index_2=double(find(symbols==3));
                index_3=double(find(symbols==2));
                symbols(index_2)= 2;
                symbols(index_3)= 3;
            end
            % 生成pam
            pamsignal = pammod(symbols,obj.TxPHY.M,0,'gray');
            pamsignal = pnorm(pamsignal);
            symbolsUp = upsample(pamsignal, obj.TxPHY.sps);
            % 生成Duob 成型器
            reverse='False';
            pulseDoub = Duob(obj.TxPHY.sps, obj.Nr.psfLength , reverse, obj.Nr.psfRollOff );
            % 生成Duob数据
            duobSignal=conv(pulseDoub,symbolsUp,'same');

        end

        % 延迟函数：
        function outSignal=delaySignal(obj,input,delay)
            switch lower(obj.Button.delaySignal)
                case 'symbol'
                    % 延迟样本数 整数 或者 小数 都能实现
                    outSignal = delay_signal(input, delay);
                case 'time'
                    % 延迟时间数
                    % delay = 50e-12;    % 50ps延迟
                    outSignal = iqdelay(input, obj.TxPHY.fs, delay);
                    outSignal=outSignal.';
            end
        end

        % 线性部分响应编码
        function [filteredPRS,prs_sig,lcoeff]=prsSignal(obj,D,taps)
            % 生成脉冲成型信号
            [~,pam_signal]=obj.dataOutput();
            %pam_signal=pnorm(pam_signal);
            % Linear encoding   1sps 选项
            lcoeff = prsFilter(D, taps, 1);
            % 归一化
            lcoeff = lcoeff / sum(lcoeff);

            % Flitering
            prs_sig = filter(lcoeff, 1, pam_signal);

            % 成型滤波器
            if strcmp(obj.Button.shapingFilter,'system')
                pulse = obj.systemHsqrt();
            elseif strcmp(obj.Button.shapingFilter,'user')
                pulse = obj.getPulse(obj.Nr.userShape);
            end
            % 脉冲成型 装载成型后的信号
            filteredPRS=obj.applyShapingFilter(prs_sig,pulse);
        end

        % 非线性部分响应编码
        function [filteredPRS,prs_sig,w]=prsSignalNonlinear(obj,D,taps1,taps2,taps3)
            % 生成脉冲成型信号
            [~,pam_signal]=obj.dataOutput();

            % Linear encoding
            lcoeff = prsFilter(D, taps1, 1);
            % 归一化
            lcoeff = lcoeff / sum(lcoeff);
            % Nonlinear encoding
            % 随机生成非线性项的权重，维度由Volterra级数的项数决定
            nlcoeff = randn(1, taps2*(taps2+1)/2 + taps3*(taps3+1)*(taps3+2)/6);
            nlcoeff = nlcoeff / sum(nlcoeff);
            % 抽头参量
            w = [lcoeff, nlcoeff];
            w=w.';
            % 一阶
            tapslen_1=taps1;
            % 二阶
            tapslen_2=taps2;
            % 三阶
            tapslen_3=taps3;
            % 前馈抽头
            x1=zeros(tapslen_1,1);

            prs_sig=zeros(length(pam_signal),1);
            sps=1; % 1sps 处理
            for i=1:length(pam_signal)
                % 按照Volterra级数形式生成非线性信号
                %一阶前馈输入
                x1 = cat(1,x1(sps+1:end),pam_signal(sps*i-sps+1:1:sps*i));
                %二阶前馈输入 % 从x1中确定好抽头的数量
                x2 = x1(round((tapslen_1-tapslen_2)/2)+1 : end - fix((tapslen_1-tapslen_2)/2));
                % 二阶核
                x2_vol = BuildVolterraInput(x2,2);
                %三阶前馈输入
                x3 = x1(round((tapslen_1-tapslen_3)/2)+1 : end - fix((tapslen_1-tapslen_3)/2));
                % 三阶核
                x3_vol =  BuildVolterraInput(x3,3);
                % 转至为行向量
                x=x1.';
                % 结合
                x_all = [x x2_vol x3_vol];
                % Flitering
                prs_sig(i) = x_all*w;
            end
            % 成型滤波器
            if strcmp(obj.Button.shapingFilter,'system')
                pulse = obj.systemHsqrt();
            elseif strcmp(obj.Button.shapingFilter,'user')
                pulse = obj.getPulse(obj.Nr.userShape);
            end
            % 脉冲成型 装载成型后的信号
            filteredPRS=obj.applyShapingFilter(prs_sig,pulse);
        end



        function [xout,coeff]=getPrsPulse(obj,taps,reverse)
            % PRS FTN格式
            %   SpS      : 每个符号的样本数（默认值=2）
            %   Nsamples : 输出脉冲的符号长度（默认值=16）
            %   reverse  : 是否翻转脉冲 ['true'/'false']（默认='false'）
            %   alpha    : 根升余弦滤波器的滚降因子（默认=0.01）
            % 输出：
            %   xout     : 双二进制脉冲成形滤波器系数
            %if nargin<2
            %    SpS=2;
            %end
            %if nargin<3
            %    Nsamples=16;
            %end
            %if nargin<4
            %    reverse='false';
            %end
            %if nargin<5
            %    alpha=0.01;
            %end
            if nargin<2
                taps=2;
            end
            if nargin<4
                reverse='false';
            end
            % 参数设置
            alpha=obj.Nr.psfRollOff;
            Nsamples=obj.Nr.psfLength;
            SpS=obj.TxPHY.sps;
            syms x;  % 定义符号变量x
            y = (1 + x)^taps;  % 构造二项式表达式 (1+x)^N
            coeff = double(flip(coeffs(expand(y))));  % 多项式展开 -> 提取系数 -> 反转顺序 -> 转为数值

            K = length(coeff) - 1;                % 多项式阶数
            max_delay = K;                         % 最大延迟符号数
            % 扩展长度（避免边界效应）
            N_extend = 2 * max_delay*SpS;   % 扩展的符号长度

            % 扩展滤波器长度（避免边界效应）
            N = Nsamples + N_extend;  % 实际生成长度 = 输出长度 + 两端各扩展SpS个样本

            % 参数说明：
            %   alpha : 滚降因子（控制带宽）
            %   N     : 滤波器长度（符号数）
            %   SpS   : 每个符号样本数
            %   'sqrt': 生成根升余弦而非普通升余弦
            p = rcosdesign(alpha,N,SpS,'sqrt');

            % 时域实现：原始脉冲 + 延迟SpS样本的脉冲
            % 向左滚动SpS 个符号

            % 初始化输出脉冲
            x = zeros(1, length(p));
            % ===== 核心：多项式加权叠加 =====
            for k = 0:K
                delay_samples = k * SpS;           % 延迟样本数
                weighted_pulse = coeff(k+1) * circshift(p, -delay_samples);
                x = x + weighted_pulse;            % 累加加权延迟脉冲
            end

            %翻转
            % 当 reverse='true' 时生成匹配滤波器
            % 脉冲翻转
            if strcmp(reverse, 'true')
                x = flip(x);
            end

            % 截取有效部分 (保留中间Nsamples长度，去除扩展的边界样本)
            start_idx = max_delay * SpS + 1;
            end_idx = length(x) - max_delay * SpS;
            xout = x(start_idx : end_idx);

        end

        % PRS Code
        function [prsSignal,coeff]=prsCodeSignal(obj,taps,reverse)
            % 生成脉冲成型信号
            [~,pam_signal]=obj.dataOutput();
            pam_signal=pnorm(pam_signal);
            [pulse,coeff]=obj.getPrsPulse(taps,reverse);
            symbolsUp = upsample(pam_signal, obj.TxPHY.sps);
            prsSignal=conv(symbolsUp,pulse,'same');
        end

        % 生成比特数
        function [data,symbols_prbs]=prbs_bits(obj)
            %参数：obj.prbsOrder，NSym，M
            %采用prbs码生成基本数据
            data = prbs1(obj.TxPHY.prbsOrder,obj.TxPHY.NSym*log2(obj.TxPHY.M),0);
            data_2bit=reshape(data,log2(obj.TxPHY.M),[]);
            %symbols_prbs = 2.^(0:log2(obj.M)-1)*data_2bit;
            symbols_prbs=double(data_2bit);
        end



        % 生成比特数
        function symbols_rand=rand_bits(obj)
            rng(obj.TxPHY.order);
            %参数：obj.prbsOrder，NSym，M
            data_2bit=randi([0,1],log2(obj.TxPHY.M),obj.TxPHY.NSym);
            % 相当于四个电平
            % bin 编码
            symbols_rand = 2.^(log2(obj.TxPHY.M)-1:-1:0)*data_2bit;
            % gray编码
            % 格雷编码映射表
            gray_map = [
                0 0;  % 二进制00 → 格雷码00 0
                0 1;  % 二进制01 → 格雷码01 1
                1 1;  % 二进制10 → 格雷码11 2
                1 0   % 二进制11 → 格雷码10 3
                ];
            % 应用格雷编码 ,2,3调换
            index_2=double(find(symbols_rand==3));
            index_3=double(find(symbols_rand==2));
            symbols_rand(index_2)= 2;
            symbols_rand(index_3)= 3;
        end

        % 异或编码
        function symbols=modcode(obj)
            %生成bit数
            switch lower(obj.Button.DataType)
                case 'prbs'
                    bits = prbs1(obj.TxPHY.prbsOrder,obj.TxPHY.NSym*log2(obj.TxPHY.M),0);
                case 'rand'
                    rng(obj.TxPHY.order);
                    %参数：obj.prbsOrder，NSym，M
                    bits=randi([0,1],log2(obj.TxPHY.M),obj.TxPHY.NSym);
            end

            c_bits=zeros(1,length(bits(:)));
            for i = 1:length(bits)
                if i==1
                    c_bits(i)=bits(i);
                else
                    %异或运算
                    c_bits(i) =xor(bits(i),c_bits(i-1));
                end
            end
            c_bits=reshape(c_bits,log2(obj.TxPHY.M),[]);
            % 二进制转换
            symbols = 2.^(log2(obj.TxPHY.M)-1:-1:0)*c_bits;
        end

        % 调制生成信号
        function pam_signal=pam(obj,symbols)
            %PAM
            % Mapeia bits to pulsos
            pam_signal = pammod(symbols,obj.TxPHY.M ,0,'gray');
        end


        %设计根升余弦脉冲成型滤波器
        function hsqrt = systemHsqrt(obj)
            hsqrt = rcosdesign(obj.Nr.psfRollOff,obj.Nr.psfLength,obj.TxPHY.sps,obj.Nr.psfShape);
        end

        % 手动设计成型滤波器
        function pulse=getPulse(obj,pulseshape)
            % pulse shape
            if strcmp(pulseshape, 'nrz')
                pulse = pulseShape('nrz', obj.TxPHY.sps);
            elseif strcmp(pulseshape, 'rrc')
                pulse = pulseShape('rrc', obj.TxPHY.sps, 4096, obj.Nr.psfRollOff, 1/obj.TxPHY.fs);
            end
            pulse = pulse / max(abs(pulse));
        end


        % 应用脉冲成型滤波器
        function    filteredSignal=applyShapingFilter(obj,symbTx,pulse)
            % Upsampling
            symbolsUp = upsample(symbTx, obj.TxPHY.sps);
            %symbolsUp = resample(symbTx, obj.TxPHY.sps,1);
            % Pulse shaping
            if strcmp(obj.Button.shapingFilter,'system')
                filteredSignal=conv(symbolsUp,pulse,'same');
            elseif strcmp(obj.Button.shapingFilter,'user')
                filteredSignal = firFilter(pulse, symbolsUp);
            end
        end


        % 对信号进行加噪
        function dataOut=addNoiseEbN0(~,datain,Eb_N0_dB)
            noise=EbN0_dB(datain,Eb_N0_dB);
            % 加入噪声
            dataOut = real(datain+noise);
        end
        % 相噪建模
        function   Pin=phaseNoise(obj,sigTx,lw)
            % 相噪建模
            phi_pn_lo = phaseNoise(lw, length(sigTx), 1/obj.TxPHY.fs);
            sigLO = exp(1i * phi_pn_lo);
            Pin=sigLO;
        end

        % 激光器建模
        function   sigLO = laserMode(obj,inputSignal,lw,RIN,Plo_dBm)
            % 转换为W
            Plo=obj.dBTow(Plo_dBm);
            % LO
            paramLO=struct();
            paramLO.P = Plo;
            paramLO.lw = lw;          % laser linewidth
            paramLO.RIN_var = RIN; % 一般设置为0
            paramLO.Fs = obj.TxPHY.fs;
            paramLO.N = length(inputSignal);
            sigLO = basicLaserModel(paramLO);
        end

        % 器件频域响应模型
        function [filt,H]=createFrequencyResponse(obj,freq,order,f3dB,verbose)
            %f3dB=20e9;
            %order=5;
            %verbose=0;
            % 归一化的截止频率
            fcnorm = f3dB/(obj.TxPHY.fs/2);
            if strcmp(obj.Implementation.responType,'Bessel')
                filt=Besslf_filter(order,fcnorm,verbose);
                % 获取频响 (已经去除延时)
                H = filt.H(freq/obj.TxPHY.fs);
            elseif strcmp(obj.Implementation.responType,'Gaussion')
                % 高斯滤波器
                nlength=length(freq);
                filt=Gaussian_filter(order,fcnorm,nlength,verbose);
                H = filt.H(freq/obj.TxPHY.fs);
            end
        end

        % ADC 加入采样时钟偏移
        function Eout=addSamplingClockOffset(obj,Ei,ppm,jitter_rms,fsBase_ADC)
            if nargin<5
                % 默认ADC的基准采样率为系统的采样率
                fsBase_ADC=obj.TxPHY.fs;
            end
            Fs_in=obj.TxPHY.fs;
            Fs_adc = fsBase_ADC*(1 + ppm/1e6);
            ppm_meas = (Fs_adc-fsBase_ADC)/(fsBase_ADC)*1e6;
            fprintf('ADC sampling rate = %.5f GS/s\n',Fs_adc/1e9);
            fprintf('ADC sampling clock drift = %.2f ppm\n',ppm_meas);

            if ~isreal(Ei)
                % Signal interpolation to the ADC's sampling frequency
                Eout = clockSamplingInterp(real(Ei), Fs_in, Fs_adc, jitter_rms) + 1i * clockSamplingInterp(imag(Ei), Fs_in, Fs_adc, jitter_rms);
            else
                % Signal interpolation to the ADC's sampling frequency
                Eout = clockSamplingInterp(Ei, Fs_in, Fs_adc, jitter_rms);
            end
        end

        % 匹配滤波
        function outSignal=matchFiltering(obj,input)
            % remove Dc
            input=input-mean(input);
            outSignal=conv(input,obj.TxPHY.hsqrt ,'same');
            outSignal=pnorm(outSignal);
        end


        % 应用信道
        function sigRxo=channelApply (obj,hch,sigTxo)

            % Channel
            % 信道响应
            % hch = [0.74 -0.514 0.37 0.216 0.062];
            % hch = [0.207, 0.815, 0.207];

            hch_up = upsample(hch, obj.TxPHY.sps);

            % 经过信道
            sigRxo=filter(hch,1,sigTxo);
            %             delay= grpdelay(hch,1,1);
            %             %去除延时点
            %             sigRxo = sigRxo(floor(delay)+1:end);

            % 使用conv进行测试
            % sigRxo=conv(sigRxo,hch,'same');
            % sigRxo=sigRxo(1:end-length(hch)+1);
            % sigRxo=sigRxo(length(hch):end);
            % 使用同步进行测试
            % [sigRxo,ref_sync,ff] = sync(sigRxo,sigRxo);
            % sigRxo=sigRxo(ff(2):end);

        end

        % 应用EA模型
        function eaOut=eaApply(obj,amptype,sigTx)
            % EA Amptype
            % amptype='Cub'; % 选用EA模型
            if strcmp(amptype,'LUT')
                % LUT的放大器模型
                pindBm =-45;                % Input power
                pin = 10.^(pindBm/10)*1e-3; %W
                paChar = pa_performance_characteristics();
                amplifier = comm.MemorylessNonlinearity('Method','Lookup table','Table',paChar,'ReferenceImpedance',50);
                sigTx = sigTx*sqrt(pin * amplifier.ReferenceImpedance);
                eaOut = real(amplifier(sigTx.'));
                plot(amplifier)
            elseif strcmp(amptype,'Cub')
                pindBm =10;      %10          % Input power
                gain = 5;        % Amplifier gain
                % ReferenceImpedance代表的是输出阻抗。。。
                amplifier = comm.MemorylessNonlinearity("Method","Cubic polynomial", ...
                    "LinearGain",gain,'IIP3',20,"AMPMConversion",10,"ReferenceImpedance",50);
                pin = 10.^(pindBm/10)*1e-3; %W
                sigTx = sigTx*sqrt(pin*amplifier.ReferenceImpedance);
                eaOut = real(amplifier(sigTx));
                plot(amplifier);
            end
            eyediagram(sigTx(1:10000),5*obj.TxPHY.sps )

        end


        % CD预补偿 FIR滤波器形式，滤波器生成
        function Pre_filter=cdPreFIR(obj,FiberLen,DER)
            % 实现了一个预均衡滤波器设计系统，通过迭代优化生成一个预加重滤波器（Pre_filter），
            % 用于在发射端预先补偿光纤传输中的色散（CD）和非线性效应。
            % 注： DER参数选取
            %DER = 1;  %% 对于GS来说需要22左右   对于FIR来说在1左右
            a=1;
            b= 10^(DER/10)*a; % 根据DER计算中间带增益
            N=900;  % 滤波器半长度

            H_FIR_Orin=sqrt([a*ones(1,N) b a*ones(1,N+1)]); % 初始频响：中间带增益√b，两侧√a
            % 迭代次数
            iteration = 40;
            C_speed = 3e8;
            % 初始响应
            H_FIR= H_FIR_Orin;
            for i=1:iteration
                % % 步骤1：发送端滤波（仅幅度影响）
                H_FIR_Tx = abs(H_FIR);
                % 步骤2：模拟光纤色散效应
                rt = FiberDispersion( H_FIR_Tx, C_speed, FiberLen, obj.TxPHY.fs   );
                % 步骤3：接收端相位响应提取
                H_FIR_Rx = abs(H_FIR_Orin).*exp(1j*angle(rt));
                % 步骤4：色散补偿计算
                H_FIR =FiberDispComp(H_FIR_Rx,C_speed, FiberLen, obj.TxPHY.fs  );
            end
            % % 去直流处理（取幅度，适用于DD系统）
            Pre_filter = abs(H_FIR).^2 - mean(abs(H_FIR).^2);
            figure;
            plot(linspace(-obj.TxPHY.fs /2,obj.TxPHY.fs /2,length(Pre_filter)),10*log10(abs(fftshift(fft(Pre_filter)))));
            title('Pre_filter');
        end

        % CD FIR Apply
        function TxWfm_PreEDC=cdFIRApply(obj,DER,FiberLen)
            % 生成
            [TxWfm,~]=obj.dataOutput();
            %DER = 1;  %% 对于GS来说需要22左右   对于FIR来说在1左右
            Pre_filter=obj.cdPreFIR(FiberLen,DER);
            % 应用
            TxWfm_PreEDC = conv(TxWfm,Pre_filter,"same");
            TxWfm_PreEDC = TxWfm_PreEDC- mean(TxWfm_PreEDC);
        end


        % CD GS
        function  TxWfm_PreEDC=cdGSApply(obj,DER,FiberLen)
            % 生成
            [TxWfm,~]=obj.dataOutput();
            % Digital extinction ratio Emulation
            % 注： DER参数选取
            %DER = 22;  %% 对于GS来说需要22左右   对于FIR来说在1左右
            TxWfm = TxWfm+ 10^(DER/20);
            % 输出数值
            disp(10*log(max(abs(TxWfm))/min(abs(TxWfm))))
            Pt=TxWfm;
            % 取幅度
            s0_t=sqrt(Pt);
            st_temp = s0_t;
            % 迭代次数
            iteration = 40;
            C_speed = 3e8;
            for i=1:iteration
                st = abs(st_temp);
                rt = FiberDispersion( st, C_speed, FiberLen, obj.TxPHY.fs  );
                rt_temp = s0_t.*exp(1j*angle(rt));
                %BW = BitRateDefault*(1+RollOff)/2+0.6e9;
                %rt_temp = Brick_wall_filter(rt_temp,BW,SampleRateDefault);
                st_temp =FiberDispComp(rt_temp,C_speed, FiberLen, obj.TxPHY.fs );
            end
            % 取功率信号
            TxWfm_PreEDC = st_temp.*conj(st_temp);
            TxWfm_PreEDC = TxWfm_PreEDC- mean(TxWfm_PreEDC);
        end


        % dBm to W (功率转换)
        function Pch = dBTow(~,Pch_dBm)
            Pch = (10 .^ (Pch_dBm / 10)) * 1e-3;
        end

        % 调制后的信号，进行功率转换
        function  signal=setSignalPower(obj,input,Channel_power_type,Pin_dBm)
            % 转换为W
            Pin=obj.dBTow(Pin_dBm);
            % 每个通道的发射功率是相对于所有偏振模式总功率，不是相对于单个偏振模式的功率。
            % 将发射功率除以偏振模式的数量可以得到每个通道在单个偏振模式下的功率。
            if strcmp(Channel_power_type,'output')
                % 这句的作用是，输出的光功率 就是 设置的功率，不用进行放大了，也即前面设置的功率为输出的光功率
                signal = sqrt(Pin) * pnorm(input);
            elseif strcmp(Channel_power_type,'input')
                % 这里的作用是， 前面设置的即为输入光功率为多少，
                signal = sqrt(Pin)* (input) ;
            end
        end

        %理想信道响应
        function  cdChannelResponse(obj,FiberLen)
            C_speed = 3e8;
            lamda = C_speed/193.1e12;
            CD_value = 17e-6 * FiberLen;
            N = 30000;
            TimeW = N/obj.TxPHY.fs;
            beta2 = -lamda^2*CD_value/2/pi/C_speed;
            w = 2*pi*[(0:N/2-1),(-N/2:-1)]/TimeW;
            % 理想的信道响应
            Response_fft=cos(-beta2*(w.^2)/2);
            % 转回时域
            Response_Tdomain=ifft(Response_fft);
            figure;
            plot(linspace(-obj.TxPHY.fs/2,obj.TxPHY.fs/2,length(Response_Tdomain)),10*log10(abs(fftshift(fft(Response_Tdomain)))));
            title('Response_Tdomain');
        end


        function [pam4demod] = pam4demod(~,sig)
            pam4demod=zeros(1,length(sig)*2);
            k = 1;
            for i=1:length(sig)
                % 'gray'编码
                switch sig(i)
                    case -3
                        pam4demod(k) =0;
                        pam4demod(k+1)=0;
                    case -1
                        pam4demod(k) = 0;
                        pam4demod(k+1) = 1;
                    case 1
                        pam4demod(k)=1;
                        pam4demod(k+1)=1;
                    case 3
                        pam4demod(k) = 1;
                        pam4demod(k+1) = 0;
                end
                k=k+2;
            end

        end

    end
end
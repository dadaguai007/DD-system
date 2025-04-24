classdef DspSyncDecoding < handle
    % 设计一个 专门 为 DSP 流程中  同步，解码 的类
    % 流程：接收信号 → 载波同步（频率补偿） → 符号同步（相位补偿） → 解调
    properties
        % 为了方便后续添加变量，采用结构体方式进行变量命名
        signalPHY; % 传入的信号参数
        Nr%无量纲参数
        Implementation;% 参考信号实施参数
    end

    methods
        function obj = DspSyncDecoding(varargin)
            %初始化类的属性
            if numel(varargin) == 6
                obj.signalPHY.fs            = varargin{1} ;% 接收信号的采样率
                obj.signalPHY.fb            = varargin{2} ;% 接收信号的波特率
                obj.signalPHY.M             = varargin{3} ;% 接收信号的格式
                obj.Nr.sps                  = varargin{4} ;% 上采样率
                obj.Nr.time_osr             = varargin{5} ;% 时钟信号的上采样率
                obj.Nr.ncut_index           = varargin{6} ;% 误码计算起始位置
                obj.Implementation.ref      = varargin{7} ;% 参考信号
            end

        end


        function clk_output=ClockRecovery(obj,input_signal)
            %数字时钟恢复（载波频率）
            f_up=obj.Nr.time_osr *obj.signalPHY.fb; % 时钟恢复所需要的采样率
            data_up = resample(input_signal,f_up,obj.signalPHY.fs);  %一般采用四倍上采样，满足后续时钟恢复频率
            % 使用时钟恢复函数
            f_clk_rec = cr(data_up.',obj.signalPHY.fb,obj.Nr.time_osr,0); % 注意：f_clk_rec是对fb的估计，不是对fb*osr的估计！
            % 创建时间轴
            [~, t_cur] = freq_time_set(length(data_up), f_up);
            %t_cur = 1/obj.Nr.time_osr/obj.signalPHY.fb*(0:length(data_up)-1); % 注意：插值前data1的采样率是osr*fb
            % 正确时间信号，重建时间轴
            t_new = 0:1/f_clk_rec/obj.Nr.time_osr:t_cur(end); % 注意：插值后data2的采样率是osr*f_clk_rec
            % 插值
            data1 = interp1(t_cur,data_up,t_new,'spline');
            % 转换到所需要的采样率
            clk_output = resample(data1,obj.Nr.sps ,obj.Nr.time_osr);
            % 归一化
            clk_output=pnorm(clk_output);
        end


        function [syncedSignal,sync_position] = Synchronization(obj, clk_output)
            % 执行同步操作
            % 参考序列 达到 需要的采样率
            ref_sps = repelem(obj.Implementation.ref,obj.Nr.sps);
            [syncedSignal,~,sync_position] = sync(clk_output,ref_sps);
        end


        function syncedSignal=Time_recovery(obj,input_signal)
            % 对输入信号执行时间（频率）恢复 和 同步
            clk_output=obj.ClockRecovery(input_signal);
            [syncedSignal,~] = obj.Synchronization(clk_output);
        end


        function   [DeWaveform,P,OptSampPhase,MaxCorrIndex]=time_phase_Recovery(obj,input_signal)
           % 参考信号
            label=obj.Implementation.ref;
            %数字时钟恢复（载波频率）
            fs_up=obj.Nr.time_osr *obj.signalPHY.fs; % 时钟恢复所需要的采样率
 [DeWaveform,P,OptSampPhase,MaxCorrIndex] = Quick_Syn_Vec(input_signal,label,1/fs_up,1/obj.signalPHY.fb  );

        end


        function [decodedData,ber] = NRZ_ExecuteDecoding(obj, eq_signal)
            % NRZ信号执行解码操作
            offset=0;
            % normalize
            sigRx_E = eq_signal - mean(eq_signal);
            sigRx_E = sigRx_E./max(sigRx_E);
            % make decision and convert into 0,1 sequence
            out = sign(sigRx_E-offset);
            decodedData = 0.5*(out+1);
            % 参考信号  重复 一定数量 ，满足解码 数量
            ref_seq=repmat(obj.Implementation.ref,1000,1);
            ref_seq=ref_seq(:);
            % NRZ 为 0,1 码元， 故解码为0,1码
            ref_seq=0.5*(ref_seq+1);

            % 解码
            [ber,num,~] = CalcBER(decodedData(obj.Nr.ncut_index:end),ref_seq(obj.Nr.ncut_index :end)); %计算误码率
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        function [decodedData,ber] = PAM_ExecuteDecoding(obj, eq_signal)
            % 针对 PAM4 解码，计算误码
            % 量化区间
            A=[-2 0 2];
            % 参考信号  重复 一定数量 ，满足解码 数量
            ref_seq=repmat(obj.Implementation.ref,1000,1);
            ref_seq=ref_seq(:);
            % 参考序列
            [~,label] = quantiz(ref_seq,A,[-3,-1,1,3]);
            label_bit=obj.pam4demod(label);
            % 接收序列
            [~,I] = quantiz(eq_signal,A,[-3,-1,1,3]);
            decodedData=obj.pam4demod(I);
            % 解码
            [ber,num,~] = CalcBER(decodedData(obj.Nr.ncut_index:end),label_bit(obj.Nr.ncut_index:end)); %计算误码率
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        function y= pam4demod(obj,sig)
            % 调制格式
            M= obj.signalPHY.M;
            % 解码类型
            y=zeros(1,length(sig)*2);
            %%PAM4
            for i = 1:length(sig)
                if sig(i) == -3
                    y(i*2-1) = 0;
                    y(i*2) = 0;
                elseif sig(i) == -1
                    y(i*2-1) = 1;
                    y(i*2) = 0;
                    %原始对应的是01
                elseif sig(i) == 1
                    y(i*2-1) = 1;
                    y(i*2) = 1;
                elseif sig(i) == 3
                    y(i*2-1) = 0;
                    y(i*2) = 1;
                    %原始对应的是10
                end
            end
            y = y(:);
        end


        function pam8demod= pam8demod(obj,sig)
            % 调制格式
            M= obj.signalPHY.M;
            % 解码类型
            pam8demod=zeros(1,length(sig)*3);
            k=1;
            %%PAM4
            for i = 1:length(sig)
                if sig(i) == -7
                    pam8demod(k) =0;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=0;
                elseif sig(i) == -5
                    pam8demod(k) =0;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=1;
                elseif sig(i) == -3
                    pam8demod(k) =0;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=1;
                elseif sig(i) == -1
                    pam8demod(k) =0;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=0;
                elseif sig(i) == 1
                    pam8demod(k) =1;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=1;
                elseif sig(i) == 3
                    pam8demod(k) =1;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=0;
                elseif sig(i) == 5
                    pam8demod(k) =1;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=0;
                elseif sig(i) == 7
                    pam8demod(k) =1;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=1;
                end
                k=k+3;
            end
        end

    end
end
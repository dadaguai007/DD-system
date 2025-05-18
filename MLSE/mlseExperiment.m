clear ;
close all;
clc;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

% make;
%% Experiments parameters
expParams = struct;
expParams.channelLengthList = 1:5;
% expParams.channelLengthList = 3;
expParams.blockLength = 16; % 全搜索只适用于示例演示！！！
expParams.expCount = 10;
expParams.snrDbRange = -10:2:10;

%% Run experiments
legendList = [];
ww = waitbar(0, "start");
% Iterate over channel lenghts
denom = length(expParams.channelLengthList)*expParams.expCount*length(expParams.snrDbRange);
fig = figure;
for channelLengthIdx = 1:length(expParams.channelLengthList)

    channelLength = expParams.channelLengthList(channelLengthIdx);
    MLSEmexFunc('new', int8(expParams.blockLength), int8(channelLength));
    % Initialize equalizer
    equalizer = mlseClasee(channelLength, expParams.blockLength);
    idx_chan = (channelLengthIdx - 1)*expParams.expCount*length(expParams.snrDbRange);
    
    % Initialize matrix for BER
    berSnrExp = zeros(length(expParams.snrDbRange), expParams.expCount);
    
    % Iterate over seed
    for expIdx = 1:expParams.expCount
        rng(expIdx);
        idx_exp = (expIdx - 1)*length(expParams.snrDbRange);
        % Rayleigh channel - CIR from CN(0,1)
        cir = 1/sqrt(2) * (randn(channelLength,1) + 1j * randn(channelLength,1));
        
        % Tx bits
        txBits = randi([0 1], expParams.blockLength, 1);
        % disp(txBits')
        
        % Modulate
        txSignal = 2*txBits-1;
        
        % Generate noise
        unitPowerNoise = 1/sqrt(2) * (randn(expParams.blockLength+channelLength-1,1) + 1j * randn(expParams.blockLength+channelLength-1,1));
        
        % Iterate over SNR
        for snrIdx = 1:length(expParams.snrDbRange)
            waitbar((idx_chan + idx_exp + snrIdx)/denom, ww, "chan " + num2str(channelLengthIdx) + " exp " + num2str(expIdx) + " snr " + num2str(expParams.snrDbRange(snrIdx)));
            snrDb = expParams.snrDbRange(snrIdx);
            
            % Received signal model
            noise = unitPowerNoise * db2mag(-snrDb);
            rxSignal = conv(txSignal, cir) + noise;
            
            % Run equalizer
            equalized = equalizer.run(rxSignal, cir);
            
            % Demodulate
            rxBits = double(equalized > 0);
            
            % Calculate BER
            berSnrExp(snrIdx, expIdx) = mean(rxBits.' ~= txBits);
        end
    end
    
    % Plot graph
    semilogy(expParams.snrDbRange, mean(berSnrExp, 2)); hold on;
    
    % Accumulate legend
    legendList = [legendList; "MLSE, CIR length = " + string(channelLength)];
end

%% Visualize results
% Calculate BPSK AWGN BER for comparison
awgnBer = berawgn(expParams.snrDbRange, 'psk', 2, 'nondiff');
% Put on plot
semilogy(expParams.snrDbRange, awgnBer);

% Accumulate legend
legendList = [legendList; "AWGN BER"];
% Put legend, apply limits, specify axis
legend(legendList); grid on; ylim([1e-5 1])
xlabel("SNR, dB"); ylabel("BER");

% Save figure
saveas(fig, "results.png");
close(ww);

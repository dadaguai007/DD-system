% PAM Amptype
amptype='Cub';

if strcmp(amptype,'LUT')
    % LUT的放大器模型
    pindBm =-45;                % Input power

    pin = 10.^(pindBm/10)*1e-3; %W
    paChar = pa_performance_characteristics();
    amplifier = comm.MemorylessNonlinearity('Method','Lookup table','Table',paChar,'ReferenceImpedance',50);
    sigTx = sigTx*sqrt(pin * amplifier.ReferenceImpedance);
    Sig_Tx = real(amplifier(sigTx.'));
    plot(amplifier)
elseif strcmp(amptype,'Cub')
    pindBm =10;      %10          % Input power
    gain = 5;        % Amplifier gain
    % ReferenceImpedance代表的是输出阻抗。。。
    amplifier = comm.MemorylessNonlinearity("Method","Cubic polynomial", ...
        "LinearGain",gain,'IIP3',20,"AMPMConversion",10,"ReferenceImpedance",50);
    pin = 10.^(pindBm/10)*1e-3; %W
    sigTx = sigTx*sqrt(pin*amplifier.ReferenceImpedance);
    Sig_Tx = real(amplifier(sigTx));
    plot(amplifier);
end
eyediagram(sigTx(1:10000),5*SpS)
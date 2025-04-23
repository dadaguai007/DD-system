% Origin version was got from Jane
% log: 2014.4.2 halfLen & quadLen are added to replace original fixed value

function [DeWaveform P OptSampPhase]=Quick_Syn(RecWaveform,label,SampleT_Osci,SampleT_AWG)

seq_len = 1e4;
SeqForSync = RecWaveform(1:seq_len);

maxPhaseError=0.01;

tSeqBefSamp=(0:numel(SeqForSync)-1)*SampleT_Osci;
tSeqBefSamp=tSeqBefSamp.';
SampPhaseVec=0:maxPhaseError:(1-maxPhaseError);
SampPhaseVec=SampPhaseVec.';
CorrResult=zeros(numel(SampPhaseVec),1);

for i=1:numel(SampPhaseVec)
    tSeqAftSamp=SampPhaseVec(i)*SampleT_AWG:SampleT_AWG:tSeqBefSamp(end);
    tSeqAftSamp=tSeqAftSamp.';
    AftSampWaveform=interp1(tSeqBefSamp,SeqForSync,tSeqAftSamp,'spline');
    CorrResult(i)=max(abs(xcorr(label,AftSampWaveform)));
end

[MaxCorr MaxCorrIndex]=max(CorrResult);

if MaxCorrIndex==1
    polyForQ=polyfit([SampPhaseVec(end);SampPhaseVec(1:2)+1],[CorrResult(end);CorrResult(1:2)],2);
%     warning('Optimum resampling phase near phase Zero');
elseif MaxCorrIndex==numel(CorrResult)
    polyForQ=polyfit([SampPhaseVec(end-1:end);SampPhaseVec(1)+1],[CorrResult(end-1:end);CorrResult(1)],2);
%     warning('Optimum resampling phase near phase One');
else
    polyForQ = polyfit(SampPhaseVec(MaxCorrIndex-1:MaxCorrIndex+1),CorrResult(MaxCorrIndex-1:MaxCorrIndex+1), 2);
end

OptSampPhase=roots(polyder(polyForQ));
OptCorr=polyval(polyForQ,OptSampPhase);
OptSampPhase=OptSampPhase-floor(OptSampPhase);
% save('OptSampPhase.mat','OptSampPhase');
tSeqBefSamp=(0:numel(RecWaveform)-1)*SampleT_Osci;
tSeqBefSamp=tSeqBefSamp.';
tSeqAftSamp=OptSampPhase*SampleT_AWG:SampleT_AWG:tSeqBefSamp(end);
tSeqAftSamp=tSeqAftSamp.';
AftSampWaveform=interp1(tSeqBefSamp,RecWaveform,tSeqAftSamp,'spline');
Corr=abs(xcorr(label,AftSampWaveform));
% plot(Corr)
peak = max(Corr);
P = find(Corr>0.95*peak);
P = flipud(numel(AftSampWaveform)-P);
P(P<1) = [];
DeWaveform=AftSampWaveform;
% close
end












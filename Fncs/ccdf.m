function [ccdfx,ccdfy] = ccdf(x,sz_window,flagPlot)
% [ccdfx,ccdfy] = ccdf(x,nSym,flagPlot)
% input:
%   x: signal waveform
%   nSym: number of symbols for PAPR calculation
%   flagPlot: flag to control if to plot the results
% output:
%   ccdfx, ccdfy: the x- and y-axis of CCDF
% copyright: Tianwai@KAIST
% version: v0.1, 2018/02/12

if nargin < 3
    flagPlot = 0;
end
ccdfx = 0:0.1:15;
papr0 = ccdfx;
N = length(x);
nn = floor(N/sz_window);
xx = x(1:nn*sz_window);
xx = reshape(xx,sz_window,nn);
for idx = 1:size(xx,2)
    papr(idx) = calc_papr(xx(:,idx));
end
for i = 1:length(papr0)
 ccdfy(i) = length(find(papr >= papr0(i)))/length(papr); % # of samples exceeding papr(i)
end
if flagPlot
   semilogy(ccdfx,ccdfy);
   xlabel('PAPR_0 (dB)');
   ylabel('Pr(PAPR > PAPR_0)');
   grid on;
end
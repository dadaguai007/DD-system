%PAM Channel
% 信道响应
hch = [0.74 -0.514 0.37 0.216 0.062];
% hch = [0.207, 0.815, 0.207];
% hch_up = upsample(hch, SpS);

sigRxo=filter(hch,1,sigRxo);
delay= grpdelay(hch,1,1);
%去除延时点
sigRxo = sigRxo(floor(delay)+1:end);

% 使用conv进行测试
% sigRxo=conv(sigRxo,hch,'same');
% sigRxo=sigRxo(1:end-length(hch)+1);
% sigRxo=sigRxo(length(hch):end);
% 使用同步进行测试
% [sigRxo,ref_sync,ff] = sync(sigRxo,sigRxo);
% sigRxo=sigRxo(ff(2):end);
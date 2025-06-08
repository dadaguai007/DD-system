% conv 与 filter 函数如此使用，效果等价
channel_result = filter(channel, 1, trSymVec);

channel_result1 = conv(trSymVec, channel);
channel_result1 = channel_result1(1:end-channel_length+1);
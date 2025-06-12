function [wout, error, w] = my_ffe_lms(win, ntaps, step, ref, sps, iter)
w = zeros(ntaps, 1);
w(round(ntaps/2)) = 1;
L = floor((length(win)-ntaps)/sps+1);
wout = zeros(L, 1);
trainlen = length(ref);
error = zeros(iter*trainlen, 1);
tempIdx = 1:ntaps;
% training
for i = 1:iter
    for j = 1:trainlen
        idx = tempIdx + (j-1)*sps;
        block = win(idx);
        temp = w.' * block;
        error((i-1)*trainlen+j) = ref(j) - temp;
        w = w + step * error((i-1)*trainlen+j) * conj(block);
    end
end
% equalizing
for i = 1:L
    idx = tempIdx + (i-1)*sps;
    block = win(idx);
    wout(i) = w.' * block;
end
end


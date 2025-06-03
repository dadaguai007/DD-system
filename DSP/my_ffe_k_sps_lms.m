function [wout, err, w] = my_ffe_k_sps_lms(win, ML, step, ref, sps, iter)
    tb = [1, 1, 1, 1, 1;
          2, 1, 1, 1, 1;
          2, 1, 2, 1, 1;
          2, 1, 2, 1, 2;
          2, 2, 2, 2, 1;
          2, 2, 2, 2, 2];
    switch sps
        case 1
            tbref = tb(1, :);
        case 1.2
            tbref = tb(2, :);
        case 1.4
            tbref = tb(3, :);
        case 1.6
            tbref = tb(4, :);
        case 1.8
            tbref = tb(5, :);
        case 2
            tbref = tb(6, :);
        otherwise
            error('Only support 1, 1.2, 1.4, 1.6, 1.8, 2sps.');
    end
    numtaps = round(ML * sps);
    w = zeros(numtaps, 5);
    w(round(numtaps/2), :) = 1;
    trainlen = length(ref);
    err = zeros(iter*trainlen, 1);
    enum = floor(length(win) / sps);
    L = enum - round(max(ML)*sps/2) + 1;
    wout = zeros(L, 1);

    block = zeros(length(w), 1);
    for i = 1:iter
        number = 1;
        for j = 1:trainlen
            index = mod(j, 5);
            if index == 0
                index = 5;
            end
            if tbref(index) == 1
                block(1:numtaps) = [win(number); block(1:numtaps-1)];
                number = number + 1;
            else
                block(1:numtaps) = [win(number+1); win(number); block(1:numtaps-2)];
                number = number + 2;
            end
            y = w(:, index).' * block;
            err((i-1)*trainlen+j) = ref(j) - y;
            w(:, index) = w(:, index) + step*err((i-1)*trainlen+j)*conj(block);
        end
    end

    block = zeros(length(w), 1);
    number = 1;
    for i = 1:L
        index = mod(i, 5);
        if index == 0
            index = 5;
        end
        if tbref(index) == 1
            block(1:numtaps) = [win(number); block(1:numtaps-1)];
            number = number + 1;
        else
            block(1:numtaps) = [win(number+1); win(number); block(1:numtaps-2)];
            number = number + 2;
        end
        wout(i) = w(:, index).' * block;
    end
end
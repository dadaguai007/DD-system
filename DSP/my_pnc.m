function wout = my_pnc(win, alpha, constellation, stage)
if nargin < 4
    stage = 2;
end
wout = win;
Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);

for i = 2:length(win)-1
    if error_judge(i-1) == error_judge(i+1)
        temp = (error(i-1)+error(i+1))/2;
    else
        temp = 0;
    end
    wout(i) = wout(i) + alpha*temp;
end

if stage >= 2
    Q = decision(wout, constellation);
    error = wout - Q;
    error_judge = sign(error);

    for i = 2:length(win)-1
        if error_judge(i-1) == error_judge(i+1)
            % temp = 0;
            temp = (error(i-1)+error(i+1));
        else
            if abs(error(i-1)) >= abs(error(i+1))
                temp = error(i+1) / 2;
            else
                temp = error(i-1) / 2;
            end
        end
        wout(i) = wout(i) + alpha*temp;
    end 
end

if stage == 3
    Q = decision(wout, constellation);
    error = wout - Q;
    error_judge = sign(error);

    for i = 2:length(win)-1
        if error_judge(i-1) == error_judge(i+1)
            temp = error(i-1)+error(i+1);
            % temp = 0;
        else
            if abs(error(i-1)) >= abs(error(i+1))
                temp = error(i+1) / 2;
            else
                temp = error(i-1) / 2;
            end
        end
        wout(i) = wout(i) + alpha*temp;
    end
end
end
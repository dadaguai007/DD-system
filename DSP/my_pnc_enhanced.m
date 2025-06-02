function wout = my_pnc_enhanced(win, alpha, constellation)

Q = decision(win, constellation);
error = win - Q;
error_judge = sign(error);
wout = win;

for i = 2:length(Q)-1
    if error_judge(i-1) == error_judge(i+1)
        if error_judge(i-1) == error_judge(i)
            temp = error(i-1)+error(i+1); % +++/--- 判决有可能发生错误，更正。
        else
            temp = 0; % +-+/-+- 证明判决正确，不更正。
        end
    else
        % temp = 0;
        judge = sign(error(i-1)+error(i+1));
        if error_judge(i) ~= judge
            temp = 0;
        else
            temp = error(i-1)+error(i+1);
        end
    end
    wout(i) = wout(i) + alpha*temp/2;
end

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 2:length(Q) - 1
    if error_judge(i-1) == error_judge(i+1)
        if error_judge(i-1) == error_judge(i)
            temp = error(i-1)+error(i+1); % +++/--- 判决有可能发生错误，更正。
        else
            temp = 0; % +-+/-+- 证明判决正确，不更正。
        end
    else
        % temp = 0;
        judge = sign(error(i-1)+error(i+1));
        if error_judge(i) ~= judge
            temp = 0;
        else
            temp = error(i-1)+error(i+1);
        end
    end
    % judge = sign(error(i-1)+error(i+1));
    % if error_judge(i) ~= judge
    %     temp = 0;
    % else
    %     temp = error(i-1)+error(i+1);
    % end
    wout(i) = wout(i) + alpha*temp/2;
end

end
classdef mlse
    % MLSE equalizer

    properties
        channelLength
        blockLength
        states_viterbi
        states_mlse
    end

    methods
        function this = mlse(channelLength, blockLength)
            this.channelLength = channelLength;
            this.blockLength = blockLength;
            this.states_mlse = this.states_calc(blockLength);
            this.states_viterbi = this.states_calc(channelLength);
        end

        function states = states_calc(varargin)
            len = varargin{2};
            states_char = dec2bin(0:2^len - 1);
            states_char_reshaped = reshape(states_char, 1, [])';
            states = reshape(str2num(states_char_reshaped), [], len);
            states = 2*states - 1;
        end

        function equalized = run(this, receivedSignal, estimatedCir)
            assert(iscolumn(receivedSignal), "Input signal must be column vector");
            assert(iscolumn(estimatedCir), "Estimated CIR must be column vector");
            assert(length(receivedSignal) == this.channelLength + this.blockLength - 1, "Input signal has wrong size");
            assert(length(estimatedCir) == this.channelLength, "Estimated CIR has wrong size");
 
            equalized = mlse_full(this, receivedSignal, estimatedCir);

            equalized_mex = MLSEmexFunc('run', this.prepare_data_for_mex(receivedSignal.'), this.prepare_data_for_mex(estimatedCir.'));
       
            if ~isequal(equalized, double(equalized_mex))
                error("MATLAB and C++ realisations show different results");
            else 
                disp("MATLAB and C++ realisations are the same");
            end 

            %             equalized = viterbi(this, receivedSignal, estimatedCir);
        end

        function data_out = prepare_data_for_mex(this, data_in)
            data_out = [real(data_in); imag(data_in)];
            data_out = reshape(data_out, 1, []);
            data_out = single(data_out);
        end

        function equalized = viterbi(this, receivedSignal, estimatedCir)
            states_modulated = (this.states_viterbi*2 - 1);
            candidates = states_modulated * estimatedCir(end:-1:1);
            metrics = zeros(size(this.states_viterbi, 1), this.blockLength + 1);
            prev_state = zeros(size(this.states_viterbi, 1),this.blockLength + 1);
            N = size(this.states_viterbi, 1);
            for i = 1:receivedSignal%this.blockLength
                if (i < this.channelLength)
                    tail_cand = states_modulated(:, end:-1:end - i + 1) * estimatedCir(1:i);
                    cur_err = abs(receivedSignal(i) - tail_cand);
                elseif(i>this.blockLength)
                    tail_cand = states_modulated(:, end:-1:end - i + 1) * estimatedCir(end:-1:i - this.blockLength + 1);
                    cur_err = abs(receivedSignal(i) - tail_cand);
                else
                    cur_err = abs(receivedSignal(i) - candidates);
                end
                prev_metr = metrics(:, i);
                prev_metr_resh = [prev_metr(1:N/2) prev_metr(N/2 + 1:end)];
                [minimal_prev_metr, idx] = min(prev_metr_resh.');
                correct_idx = (idx - 1) * N/2 + (1:N/2);
                prev_state(:, i) = reshape([correct_idx;correct_idx], N, []);
                tmp = reshape([minimal_prev_metr; minimal_prev_metr], N, []);
                metrics(:, i + 1) = tmp + cur_err;
            end
            equalized = zeros(1, this.blockLength);
            [~, idx] = min(metrics(:, end));
            equalized(end) = mod(idx+1, 2);
            for i = this.blockLength - 1:-1:2
                idx = prev_state(idx, i);
                equalized(i) =  mod(idx+1,2);
            end
            a = 1;
        end

        function equalized = mlse_full(this, receivedSignal, estimatedCir)
            candidates = zeros(size(this.states_mlse, 1), this.channelLength + this.blockLength - 1);
            for i = 1:size(this.states_mlse, 1)
                candidates(i, :) = conv(this.states_mlse(i, :), estimatedCir);
            end
            err = sum(abs(candidates - receivedSignal.'), 2);
            [~, idx] = min(err);
            equalized = this.states_mlse(idx, :);
        end
    end
end


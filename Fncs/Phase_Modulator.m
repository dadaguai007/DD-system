function Ao = Phase_Modulator(Ai, u, Vb,Vpi)
% Optical Phase Modulator (PM).

% Ensure u is a column vector
if size(u, 1) > 1
    u = u';
end

% Match dimensions of Ai and u
if isscalar(Ai)
    Ai = Ai * ones(size(u));
end


% Calculate modulated optical field
Ao = Ai .* exp(1i * (u / Vpi) * pi+1i*Vb*Vpi/2);
end

function [ Out ] = FiberDispersion( In, C_speed, FiberLen, SampleRateDefault  )

% C_speed = 3e8;
% FiberLen = 0*80e3;
% SampleRateDefault = BaudRate*Up;
lamda = C_speed/193.1e12;
CD_value = 17e-6 * FiberLen;
N = length(In);
TimeW = N/SampleRateDefault;
beta2 = -lamda^2*CD_value/2/pi/C_speed;
w = 2*pi*[(0:N/2-1),(-N/2:-1)]/TimeW;
InFFT = fft(In);
OutFFT = InFFT.*exp(1i*-beta2*(w.^2)/2);
Out = ifft(OutFFT);

end


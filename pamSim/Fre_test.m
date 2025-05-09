% 频域滤波的大致处理方式
for i=1:L
input_x=x;
x_hat=buffer(input_x);
X=fft(x_hat);
X_hat=conj(X);
Y=X*w;


y=ifft(Y);
e=d-y;
e_hat=[zeros(1,N-1),e];
E=fft(e_hat);

T=X_hat*E;
w=w+T*2*miu;
end
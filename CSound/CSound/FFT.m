% Returns the single sided FFT and the freq range To be used
% complementarily with FR2TR

function [f,Y]=FFT(T,x,string)

yfft=fftshift(fft(x));



fftyscaling=length(T)/2;
fftxaxis=linspace(-1/2/(T(2)-T(1)),1/2/(T(2)-T(1)),length(T));
fftxaxis_ss=fftxaxis(length(fftxaxis)/2+1:end);

yfft_ss=yfft(floor(length(yfft)/2)+1:end)./fftyscaling;

if(strcmp(string,'display')==1)
figure
plot(fftxaxis_ss,abs(yfft_ss));

figure
plot(fftxaxis_ss,unwrap(phase(yfft_ss))./2/pi);
end
f=fftxaxis_ss;
Y=conj(yfft_ss);                % exp(-iwt) convention
end

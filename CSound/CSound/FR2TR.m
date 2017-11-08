% This function takes in a Freq vector and the TF of  a Linear System and
% constructs the Time domain response.
% Takes in FFT of exp(-i omega t)


function [t_axis,signal]=FR2TR(Freq,H,F_stim,disp_flag)

t_axis=linspace(0,1./(Freq(2)-Freq(1)),length(Freq));


%--------------------------Stimulus frequency spectrum------------
%
freq_indx=find(Freq>F_stim,1);

if(F_stim~=0)
    F_signal=max(Freq)*(Freq==Freq(freq_indx));  % For delta function                              % Frequency spectrum of the input
    
else
    F_signal=ones(1,length(Freq)).';                                                         % For impulse response
%     fprintf('\nNo stimulation!\n');
end
%-----------------------------------------------------------------

%--------------------------Response spectrum----------------------

FRF=H.*F_signal;                                                    % Frequency spectrum of the output

%------------------------------------------------------------------


signal = length(Freq)*ifft(conj(FRF)./2,'symmetric');    % assume exp(-iwt) convention

if (strcmp(disp_flag,'display')==1)
figure;plot(t_axis,signal);
end
end



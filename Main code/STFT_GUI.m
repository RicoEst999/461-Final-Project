function stft_gui_slider()
%loads audio file and displaays:
% - spectogram (with gui that adjust the window by 32 or just a slider)
% - Time domain waveform of the original and cleaned signal 
% - applying a magnitude threshold masking to STFT
%
%

%%note to use the function outside there will be a need to use paraemters
%% like function [x, fs, clean_audio] = stft_gui_slider() and declare outsid ethe fucntion - clean_audio = [];


    [x, fs] = audioread('sent001.wav'); 
    x = x(:,1);  
    clean_audio = [];
    default_win = round(.02 * fs);  %making it default window of 20ms
    min = 64;
    max = 1024;

    %Figure 1: Spectrogram
    fig1 = figure('Name','STFT Window Slider','Position',[200, 200, 800, 600]);
    axSpec = axes('Parent', fig1);

    %slider that moves from 64 to 1024, starts at the default window
    uicontrol('Parent', fig1, 'Style','text','String','Window Size:', ...
              'Units','normalized','Position',[.05 .02 .15 .05]);
    slider = uicontrol('Parent', fig1, 'Style','slider', ...
              'Min', min, 'Max', max, 'Value', default_win, ...
              'SliderStep', [32/(max-min) (2*32)/(max-min)], ...
              'Units','normalized', 'Position',[0.32 0.02 0.6 0.05], ...
              'Callback', @sliderChanged);

    % Figure 2: Time-Domain Plots
    fig2 = figure('Name','Time Domain Audio','Position',[1020, 200, 800, 600]);
    Orig = subplot(2,1,1, 'Parent', fig2);
    Clean = subplot(2,1,2, 'Parent', fig2);

    % original waveform plot
    t = (0:length(x)-1)/fs;
    plot(Orig, t, x);
    title(Orig, 'Original Audio Waveform');
    xlabel(Orig, 'Time (s)');
    ylabel('Amplitude');

    %%%%
    updateSTFT(round(slider.Value), false);

    % Slider callback (with sound)
    function sliderChanged(~, ~)
        win = round(slider.Value);
        updateSTFT(win, true);  % true = play sound
    end

    % Core plotfunction
    function updateSTFT(win, ~)
        overlap = round(win / 2);
        nfft = 1024;

        [S,F,T] = spectrogram(x, hamming(win), overlap, nfft, fs);
        mag = abs(S);%computing the magnitude of STFT
        %based mask
        threshold = median(mag(:)) * 0.5;
        mask = mag > threshold;
        masked_S = S .* mask;

        % Plot spectrogram (Figure 1)
        figure(fig1); axes(axSpec); cla;
        imagesc(T, F, 20*log10(abs(S))); 
        axis xy;
        title(sprintf('Spectrogram (Window = %d samples)', win));
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;

        % inverse of stft - ISTFT
        clean_audio = istft(masked_S, fs, ...
            'Window', hamming(win), ...
            'OverlapLength', overlap, ...
            'FFTLength', nfft);
        clean_audio = real(double(clean_audio));

        % Plot cleaned audio (Figure 2)
        t_clean = (0:length(clean_audio)-1)/fs;
        figure(fig2); axes(Clean); cla;
        plot(t_clean, clean_audio);
        title('Cleaned Audio After Masking');
        xlabel('Time (s)');
        ylabel('Amplitude');
         
    end
end


% [x, fs, cleaned] = stft_gui_slider();
% disp('Playing original audio...');
% sound(x, fs);
% pause(length(x)/fs + 0.5);
% 
% disp('Playing cleaned audio...');
% sound(cleaned, fs);
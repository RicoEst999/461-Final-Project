function Combined_Filter_GUI()
% BANDPASS_GUI_FILTER
% Interactive GUI for applying a band-pass filter to a speech signal.
% Lets the user adjust low and high cutoff frequencies via sliders
% and hear the effect in real time.
% 
% VARIABLE LIST BY MODULE:
%
% Main Function: Combined_Filter_GUI()
%   filtered_audio : Audio signal after bandpass filtering (output from bandpass_gui)
%   fs             : Sampling frequency of the audio file
%
% Bandpass GUI Function: bandpass_gui()
%   x                : Original input audio signal loaded from file
%   fs               : Sampling frequency of the audio signal
%   t                : Time axis vector for plotting the waveform
%   f                : Handle to the GUI figure for the bandpass filter
%   ax               : Handle to the axes for plotting the waveform
%   h_plot           : Handle to the plot object, used for updating the waveform display
%   default_low_cut  : Default low cutoff frequency for the bandpass filter (Hz)
%   default_high_cut : Default high cutoff frequency for the bandpass filter (Hz)
%   s1               : Handle to the slider for adjusting the low cutoff frequency
%   s2               : Handle to the slider for adjusting the high cutoff frequency
%   t1               : Handle to the text label displaying the low cutoff frequency value
%   t2               : Handle to the text label displaying the high cutoff frequency value
%   b                : Numerator coefficients of the Butterworth bandpass filter
%   a                : Denominator coefficients of the Butterworth bandpass filter
%   filtered_audio   : Audio signal after applying the bandpass filter (output of this function)
%   low_cut          : Current low cutoff frequency from slider s1 (within updateFilter)
%   high_cut         : Current high cutoff frequency from slider s2 (within updateFilter)
%   Wn               : Normalized cutoff frequencies for butter() (within updateFilter)
%
% STFT GUI Function: stft_gui_slider(x, fs)
%   x              : Input audio signal (this is `filtered_audio` from bandpass_gui)
%   fs             : Sampling frequency of the input audio signal
%   clean_audio    : Denoised audio signal after STFT processing and ISTFT
%   default_win    : Default window size for STFT in samples (approx. 30 ms)
%   min_win        : Minimum allowed window size for STFT (samples)
%   max_win        : Maximum allowed window size for STFT (samples, clamped by signal length or 1024)
%   fig1           : Handle to the GUI figure for displaying the spectrogram
%   axSpec         : Handle to the axes for the spectrogram plot
%   slider         : Handle to the UI control (slider) for adjusting the STFT window size
%   fig2           : Handle to the GUI figure for displaying time-domain audio plots
%   Orig           : Handle to the subplot for the original (input to STFT) waveform
%   Clean          : Handle to the subplot for the cleaned (output of STFT) waveform
%   t              : Time axis vector for the original audio waveform plot
%   win            : Current STFT window size in samples, from slider (within updateSTFT)
%   maxAllowedWin  : Maximum allowed window size to prevent errors (within updateSTFT)
%   nfft           : FFT length for the spectrogram and ISTFT (next power of 2)
%   overlap        : Hop size for STFT/ISTFT, typically 50% of the window size
%   S              : STFT matrix of the input signal
%   F              : Frequency vector for the spectrogram
%   T              : Time vector for the spectrogram
%   mag            : Magnitude of the STFT matrix (abs(S))
%   threshold      : Threshold value for magnitude masking (based on median of magnitudes)
%   mask           : Binary mask created by comparing magnitudes to the threshold
%   masked_S       : STFT matrix after applying the binary mask
%   t_clean        : Time axis vector for the cleaned audio waveform plot

    % PART 1: BANDPASS GUI 
    [filtered_audio, fs] = bandpass_gui();  % Run and get output

    %debugging 
if isempty(filtered_audio) || ~isvector(filtered_audio)
        error('Filtered audio is empty or invalid.');
end

    disp("DEBUG: filtered_audio length = " + num2str(length(filtered_audio)));
    stft_gui_slider(filtered_audio, fs);  % Launch STFT GUI
end

%% --------- BANDPASS GUI FUNCTION -----------
function [filtered_audio, fs] = bandpass_gui()
    %% combined audio 
<<<<<<< HEAD:Main code/Combined_Filter_GUI.m
    [x, fs] = audioread('Police Siren.wav');%loading first audio 
=======
    [x, fs] = audioread('white.wav');%loading first audio 
>>>>>>> 610a1fd1a03bc6ab575e2a7e5d0e62e1fad0da2f:Main code/Combined_Filter_GUI.txt
    x = x(:, 1);% Use only one channel if stereo (mono)

    t = (0:length(x)-1) / fs;% Time axis for plotting
    % Set default cutoff frequencies (typical speech band)
    default_low_cut = 300;
    default_high_cut = 3400;

       %% Create GUI Figure
    f = figure('Name','Band-Pass Filter GUI', ...
               'Units','normalized', ...
               'Position',[0.2 0.2 0.6 0.6], ...
               'Resize','on');

    %% Axes for plotting waveform
    ax = axes('Parent', f, 'Units','normalized', 'Position',[0.1 0.45 0.8 0.45]);
    h_plot = plot(ax, t, x);              % Plot the original waveform
    title(ax, 'Filtered Signal');
    xlabel(ax, 'Time (seconds)');
    ylabel(ax, 'Amplitude (normalized)');

    %% --- Low Cutoff Frequency Slider and Label ---
    uicontrol(f, 'Style','text','Units','normalized', ...
              'Position',[0.05 0.36 0.2 0.035], 'String','Low Cutoff (Hz)');
    s1 = uicontrol(f, 'Style','slider','Units','normalized', ...
                   'Min', 100, 'Max', fs/2 - 200, 'Value', default_low_cut, ...
                   'Position',[0.05 0.32 0.35 0.04], ...
                   'Callback',@updateFilter);  % Call updateFilter when moved
    t1 = uicontrol(f, 'Style','text','Units','normalized', ...
                   'Position',[0.41 0.32 0.1 0.04], ...
                   'String', sprintf('%.0f Hz', default_low_cut));

    %% --- High Cutoff Frequency Slider and Label ---
    uicontrol(f, 'Style','text','Units','normalized', ...
              'Position',[0.55 0.36 0.2 0.035], 'String','High Cutoff (Hz)');
    s2 = uicontrol(f, 'Style','slider','Units','normalized', ...
                   'Min', 300, 'Max', fs/2 - 1, 'Value', default_high_cut, ...
                   'Position',[0.55 0.32 0.35 0.04], ...
                   'Callback',@updateFilter);  % Call updateFilter when moved
    t2 = uicontrol(f, 'Style','text','Units','normalized', ...
                   'Position',[0.91 0.32 0.07 0.04], ...
                   'String', sprintf('%.0f Hz', default_high_cut));

    %% --- Play Filtered Audio Button ---
    uicontrol(f, 'Style','pushbutton','Units','normalized', ...
              'String','Play Filtered Audio', ...
              'Position',[0.35 0.22 0.3 0.05], ...
              'Callback',@playFiltered);

    %% --- Reset Button to Restore Default Cutoff Frequencies ---
    uicontrol(f, 'Style','pushbutton','Units','normalized', ...
              'String','Reset to Default', ...
              'Position',[0.35 0.15 0.3 0.05], ...
              'Callback',@resetDefaults);

    %% --- Static Info Text on Human Speech Range ---
    uicontrol(f, 'Style','text','Units','normalized', ...
              'Position',[0.2 0.06 0.6 0.04], ...
              'FontSize', 10, ...
              'ForegroundColor', [0 0 0.5], ...
              'String','Typical human speech frequency range: 300 Hz – 3400 Hz');
    
    %% Initialize output signal
    filtered_audio = x;

    %% --- CALLBACK: Update the Filter and Plot ---
    function updateFilter(~, ~)
        % Get slider values
        low_cut = s1.Value;
        high_cut = s2.Value;

        % Update text labels with slider values
        t1.String = sprintf('%.0f Hz', low_cut);
        t2.String = sprintf('%.0f Hz', high_cut);

        % Ensure at least 50 Hz gap between low and high cutoff
        if high_cut <= low_cut + 50
            high_cut = low_cut + 50;
            s2.Value = high_cut;
            t2.String = sprintf('%.0f Hz', high_cut);
        end

        % Normalize cutoff frequencies for Butterworth filter design
        Wn = [low_cut high_cut] / (fs/2);
        Wn = max(min(Wn, 0.999), 0.001);% Clamp values to avoid invalid cutoff
        try
            % Design a 4th-order Butterworth bandpass filter
            [b, a] = butter(4, Wn, 'bandpass');
            filtered_audio = filter(b, a, x); % Apply filter to signal
            set(h_plot, 'YData', filtered_audio); % Update plot
            drawnow;
        catch err
            disp('Error designing filter:');
            disp(err.message);
        end
    end

%% --- CALLBACK: Play the Filtered Audio ---
    function playFiltered(~, ~)
        sound(filtered_audio, fs);
    end

 %% --- CALLBACK: Reset to Default Frequency Settings ---
    function resetDefaults(~, ~)
        s1.Value = default_low_cut;
        s2.Value = default_high_cut;
        t1.String = sprintf('%.0f Hz', default_low_cut);
        t2.String = sprintf('%.0f Hz', default_high_cut);
        updateFilter();
    end
updateFilter();  % ensure initial filter is applied
    uiwait(f);  % Pause until GUI is closed
end

%% --------- STFT GUI FUNCTION -----------
function stft_gui_slider(x, fs)
    % loads audio (x,fs) and displays:
    %  - spectrogram (with a slider over window size)
    %  - original and “cleaned” time‐domain waveforms
    %  - uses a 0.5×median‐based soft/Wiener‐type mask so speech isn't choppy

    default_win = 1024;                 % default window length in samples
    min_win     = 256;
    max_win     = min(2048, length(x));

    % FIGURE 1: Spectrogram + Slider
    fig1   = figure('Name','STFT Window Slider','Position',[200,200,800,600]);
    axSpec = axes('Parent', fig1);

    % Slider label
    uicontrol('Parent', fig1, 'Style','text', ...
              'String','Window Size (samples):', ...
              'Units','normalized', 'Position',[0.05 0.02 0.20 0.05]);

    % Slider (continuous, but will round to nearest multiple of 32 below)
    slider = uicontrol('Parent', fig1, 'Style','slider', ...
              'Min', min_win, 'Max', max_win, 'Value', default_win, ...
              'SliderStep',[32/(max_win-min_win)  (32)/(max_win-min_win)], ...
              'Units','normalized', 'Position',[0.32 0.02 0.60 0.05], ...
              'Callback', @sliderChanged);

    % FIGURE 2: Time‐Domain Plots
    fig2 = figure('Name','Time Domain Audio','Position',[1020,200,800,600]);
    Orig = subplot(2,1,1, 'Parent', fig2);
    Clean = subplot(2,1,2, 'Parent', fig2);

    % Plot original waveform once
    t = (0:length(x)-1)/fs;
    plot(Orig, t, x);
    title(Orig,'Original Audio Waveform');
    xlabel(Orig,'Time (s)');
    ylabel(Orig,'Amplitude');

    updateSTFT();

<<<<<<< HEAD:Main code/Combined_Filter_GUI.m
    % Slider callback (with sound)
    function sliderChanged(~, ~)
        updateSTFT(round(slider.Value), false);
    end

    % Core plotfunction
    function updateSTFT(varargin)
    win = round(slider.Value);
    maxAllowedWin = 4096;        % <-- change as needed
    win = min(win, maxAllowedWin);

       
       nfft = 2^( nextpow2( 2*(win-1) ) );  % Ensure FFT length is valid
        overlap = round(win / 2);

            [S,F,T] = spectrogram(x, hamming(win), overlap, nfft, fs);
     

        mag = abs(S);%computing the magnitude of STFT
        %based mask
        threshold =0.8* median(mag(:));
        mask = mag > threshold;
        masked_S = S .* mask;

        % Plot spectrogram (Figure 1)
        figure(fig1); 
        axes(axSpec); cla(axSpec);
        imagesc(T, F, 20*log10(abs(S))); 
        axis xy;
        title(sprintf('Spectrogram (Window = %d samples)', win));
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;
         
       clean_audio = istft(masked_S, fs, ...
    'Window',            hamming(win),  ...
    'OverlapLength',     overlap,       ...
    'FFTLength',         nfft,          ...
    'ConjugateSymmetric', true);

        clean_audio = real(double(clean_audio));
       
        % Plot cleaned audio (Figure 2)
        t_clean = (0:length(clean_audio)-1)/fs;
        figure(fig2); 
        axes(Clean); cla;
        plot(t_clean, clean_audio);
        title('Cleaned Audio After Masking');
        xlabel('Time (s)');
        ylabel('Amplitude');
         
        sound(clean_audio, fs);
     
=======
    % Slider callback simply calls updateSTFT
    function sliderChanged(~,~)
        updateSTFT();
    end

    function updateSTFT()
        % slider and snap to multiple of 32 samples
        raw_win = round(slider.Value);
        win = 32 * round(raw_win/32);
        win = max(min(win, max_win), min_win);

        % nfft 
        nfft = 2^nextpow2(2*(win - 1));

        % 75% overlap
        overlap = round(0.75 * win);

        % one‐sided STFT on x
        [S, F, T] = spectrogram(x, hamming(win), overlap, nfft, fs);

        % Soft/Wiener‐type mask (0.5×median)
        mag = abs(S);
        thr = 0.5 * median(mag(:));
        gain = (mag.^2) ./ (mag.^2 + thr^2 + eps);  
        S_mask = S .* gain;

        % spectrogram plot (magnitude in dB)
        figure(fig1);
        axes(axSpec); cla(axSpec);
        imagesc(T, F, 20*log10(mag + eps));
        axis(axSpec,'xy');
        title(axSpec, sprintf('Spectrogram (Window = %d samples)', win));
        xlabel(axSpec,'Time (s)');
        ylabel(axSpec,'Frequency (Hz)');
        colorbar(axSpec);

        % Inverse STFT (one‐sided) to get clean_audio
        clean_audio = istft(S_mask, fs, ...
              'Window',            hamming(win), ...
              'OverlapLength',     overlap, ...
              'FFTLength',         nfft, ...
              'ConjugateSymmetric', true);
        clean_audio = real(double(clean_audio));

        % cleaned waveform plot
        t_clean = (0:length(clean_audio)-1)/fs;
        figure(fig2);
        axes(Clean); cla(Clean);
        plot(Clean, t_clean, clean_audio);
        title(Clean,'Cleaned Audio After Soft Masking');
        xlabel(Clean,'Time (s)');
        ylabel(Clean,'Amplitude');

        %Playback
        sound(clean_audio, fs);
>>>>>>> 610a1fd1a03bc6ab575e2a7e5d0e62e1fad0da2f:Main code/Combined_Filter_GUI.txt
    end
end




<<<<<<< HEAD:Main code/Combined_Filter_GUI.m
%okay so now window error isn't going to happen, we are now able to get all
%STFT options, for a clear audio we want the spectrogramt o be all blue but
%again i don't think we can create a perfect clear audio and mostly with the
% motor and the low men voice just too much of the same frequency 
=======
 
>>>>>>> 610a1fd1a03bc6ab575e2a7e5d0e62e1fad0da2f:Main code/Combined_Filter_GUI.txt

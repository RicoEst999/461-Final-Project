function Combined_Filter_GUI()
%==============================
%Main Function: Combined_Filter_GUI
%==============================
%filtered_audio  : The output signal after band-pass filtering; passed to STFT GUI.
%fs              : Sampling frequency of the audio file (in Hz).
%
%==============================
%Bandpass GUI Function: bandpass_gui()
%==============================
%x               : Original input audio signal loaded from 'white_noise.wav'.
%fs              : Sampling frequency of the input signal.
%t               : Time vector for plotting the waveform.
%default_low_cut : Default low cutoff frequency for band-pass filter (300 Hz).
%default_high_cut: Default high cutoff frequency for band-pass filter (3400 Hz).
%f               : Handle for the band-pass GUI figure.
%ax              : Handle for the waveform axes in the GUI.
%h_plot          : Plot handle for updating waveform in GUI.
%s1              : Slider UI control for adjusting low cutoff frequency.
%t1              : Text label displaying value of s1 (low cutoff).
%s2              : Slider UI control for adjusting high cutoff frequency.
%t2              : Text label displaying value of s2 (high cutoff).
%low_cut         : Current low cutoff frequency selected via s1.
%high_cut        : Current high cutoff frequency selected via s2.
%Wn              : Normalized cutoff frequencies for Butterworth filter.
%b, a            : Filter coefficients for 4th-order Butterworth filter.
%
%==============================
%STFT GUI Function: stft_gui_slider(x, fs)
%==============================
%x               : Audio input to STFT GUI (filtered signal).
%fs              : Sampling frequency of the signal.
%clean_audio     : Output signal after STFT masking and ISTFT reconstruction.
%default_win     : Default STFT window size in samples (~30 ms).
%min_win         : Minimum allowed STFT window size (64 samples).
%max_win         : Maximum allowed STFT window size (≤ 1024 or signal length).
%fig1            : Figure handle for spectrogram plot.
%axSpec          : Axes handle for displaying the spectrogram.
%slider          : Slider UI control for adjusting STFT window size.
%fig2            : Figure handle for time-domain waveform plots.
%Orig            : Subplot handle for original waveform plot.
%Clean           : Subplot handle for cleaned waveform plot.
%t               : Time vector for plotting the original waveform.
%win             : Current STFT window size (updated via slider).
%overlap         : Overlap between adjacent STFT frames (50% of window).
%nfft            : FFT length for STFT and ISTFT.
%S               : STFT matrix of the input signal.
%F               : Frequency vector for spectrogram.
%T               : Time vector for spectrogram.
%mag             : Magnitude of STFT matrix (abs(S)).
%threshold       : Threshold for binary masking (median-based).
%mask            : Binary mask identifying significant STFT bins.
%masked_S        : STFT matrix after applying binary mask.
%t_clean         : Time vector for plotting cleaned signal waveform.
%doPlay          : Boolean flag to control whether audio plays automatically.


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
    [x, fs] = audioread('white_noise.wav');% Load noisy speech file
    x = x(:, 1);% Use only one channel if stereo
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
%loads audio file and displaays:
% - spectogram (with gui that adjust the window by 32 or just a slider)
% - Time domain waveform of the original and cleaned signal 
% - applying a magnitude threshold masking to STFT

    clean_audio = [];
    default_win = round(.03 * fs);  %making it default window of 30ms
    min_win = 64;
    max_win = min(1024, length(x));

    %Figure 1: Spectrogram
    fig1 = figure('Name','STFT Window Slider','Position',[200, 200, 800, 600]);
    axSpec = axes('Parent', fig1);

    %slider that moves from 64 to 1024, starts at the default window
    uicontrol('Parent', fig1, 'Style','text','String','Window Size:', ...
              'Units','normalized','Position',[.05 .02 .15 .05]);
    slider = uicontrol('Parent', fig1, 'Style','slider', ...
              'Min', min_win, 'Max', max_win, 'Value', default_win, ...
              'SliderStep', [32/(max_win-min_win) (2*32)/(max_win-min_win)], ...
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

    updateSTFT(round(slider.Value), false);

    % Slider callback (with sound)
    function sliderChanged(~, ~)
        updateSTFT(round(slider.Value), true);
    end

    % Core plotfunction
    function updateSTFT(win, doPlay)
            if isempty(x) || length(x) < 64
        warning('Signal too short.');
        return;
            end
            win = round(win);
    win = min(win, length(x));

        overlap = round(win / 2);
       nfft = max(1024, 2^nextpow2(win));   % Ensure FFT length is valid

       try
            [S,F,T] = spectrogram(x, hamming(win), overlap, nfft, fs);
        catch err
            warning("Spectrogram failed: " + err.message);
            return;
       end 

        if size(S,1) < win
        warning('STFT result too short for ISTFT.');
        return;
        end

        mag = abs(S);%computing the magnitude of STFT
        %based mask
        threshold = median(mag(:)) * 0.5;
        mask = mag > threshold;
        masked_S = S .* mask;

        % Plot spectrogram (Figure 1)
        figure(fig1); 
        axes(axSpec); cla;
        imagesc(T, F, 20*log10(abs(S))); 
        axis xy;
        title(sprintf('Spectrogram (Window = %d samples)', win));
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;

          try
        clean_audio = istft(masked_S, fs, ...
            'Window', hamming(win), ...
            'OverlapLength', overlap, ...
            'FFTLength', nfft);
        clean_audio = real(double(clean_audio));
    catch err
        warning("ISTFT failed: " + err.message);
        return;
          end
       
        % Plot cleaned audio (Figure 2)
        t_clean = (0:length(clean_audio)-1)/fs;
        figure(fig2); 
        axes(Clean); cla;
        plot(t_clean, clean_audio);
        title('Cleaned Audio After Masking');
        xlabel('Time (s)');
        ylabel('Amplitude');
         
     if doPlay
        sound(clean_audio, fs);
     end
    end
end


%when the error appear instead of crashing we get a warnning that the STFT
%dont have enough freq bins to construct a signal via ISTFT, then we we can
%try other values to safely choose other options for window size
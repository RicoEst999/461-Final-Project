function Bandpass_GUI_Filter
    % BANDPASS_GUI_FILTER
    % Interactive GUI for applying a band-pass filter to a speech signal.
    % Lets the user adjust low and high cutoff frequencies via sliders
    % and hear the effect in real time.

    % === Variable Key ===
    % x, fs       : Audio signal and sampling rate
    % t           : Time vector for plotting
    % s1, s2      : Handles for low and high cutoff sliders
    % t1, t2      : Text boxes showing selected cutoff frequencies
    % ax          : Axes object for waveform display
    % h_plot      : Handle for plot of the filtered waveform
    % filtered_audio : The signal after filtering
    % default_low_cut, default_high_cut : Initial cutoff values (speech band)

    %% Load Audio File
    [x, fs] = audioread('noisy_speech.wav');  % Load noisy speech file
    x = x(:, 1);                              % Use only one channel if stereo
    t = (0:length(x)-1) / fs;                 % Time axis for plotting

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
            if high_cut >= fs/2 - 1
                high_cut = fs/2 - 1;
            end
            s2.Value = high_cut;
            t2.String = sprintf('%.0f Hz', high_cut);
        end

        % Normalize cutoff frequencies for Butterworth filter design
        Wn = [low_cut high_cut] / (fs/2);
        Wn = max(min(Wn, 0.999), 0.001);  % Clamp values to avoid invalid cutoff

        try
            % Design a 4th-order Butterworth bandpass filter
            [b, a] = butter(4, Wn, 'bandpass');
            filtered_audio = filter(b, a, x);      % Apply filter to signal
            set(h_plot, 'YData', filtered_audio);  % Update plot
            drawnow;
        catch err
            disp('Error designing filter:');
            disp(err.message);
        end
    end

    %% --- CALLBACK: Play the Filtered Audio ---
    function playFiltered(~, ~)
        sound(filtered_audio, fs);  % Play current filtered signal
    end

    %% --- CALLBACK: Reset to Default Frequency Settings ---
    function resetDefaults(~, ~)
        s1.Value = default_low_cut;  % Reset low cutoff slider
        s2.Value = default_high_cut; % Reset high cutoff slider
        t1.String = sprintf('%.0f Hz', default_low_cut);
        t2.String = sprintf('%.0f Hz', default_high_cut);
        updateFilter();              % Re-apply filter with default settings
    end
end

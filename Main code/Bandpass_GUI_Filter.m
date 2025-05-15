function Bandpass_GUI_Filter
    % === Variable Key ===
    % x, fs: Signal and sampling rate
    % t: Time vector
    % s1, s2: Slider handles
    % t1, t2: Text labels for frequencies
    % ax: Plot axes handle
    % h_plot: Plot object
    % filtered_audio: Filtered output signal
    % default_low_cut, default_high_cut: Starting slider values

    %% Load Audio File
    [x, fs] = audioread('noisy_speech.wav');
    x = x(:, 1);
    t = (0:length(x)-1) / fs;

    % Default filter settings
    default_low_cut = 300;
    default_high_cut = 3400;

    % Create figure
    f = figure('Name','Band-Pass Filter GUI',...
               'Units','normalized',...
               'Position',[0.2 0.2 0.6 0.6],...
               'Resize','on');

    % Axes for waveform plot
    ax = axes('Parent', f, 'Units','normalized', 'Position',[0.1 0.45 0.8 0.45]);
    h_plot = plot(ax, t, x);
    title(ax, 'Filtered Signal');
    xlabel(ax, 'Time (seconds)');
    ylabel(ax, 'Amplitude (normalized)');

    % --- Low Cutoff Slider and Label ---
    uicontrol(f, 'Style','text','Units','normalized',...
              'Position',[0.05 0.36 0.2 0.035],'String','Low Cutoff (Hz)');
    s1 = uicontrol(f, 'Style','slider','Units','normalized',...
                   'Min', 100, 'Max', fs/2 - 200, 'Value', default_low_cut,...
                   'Position',[0.05 0.32 0.35 0.04],...
                   'Callback',@updateFilter);
    t1 = uicontrol(f, 'Style','text','Units','normalized',...
                   'Position',[0.41 0.32 0.1 0.04],...
                   'String', sprintf('%.0f Hz', default_low_cut));

    % --- High Cutoff Slider and Label ---
    uicontrol(f, 'Style','text','Units','normalized',...
              'Position',[0.55 0.36 0.2 0.035],'String','High Cutoff (Hz)');
    s2 = uicontrol(f, 'Style','slider','Units','normalized',...
                   'Min', 300, 'Max', fs/2 - 1, 'Value', default_high_cut,...
                   'Position',[0.55 0.32 0.35 0.04],...
                   'Callback',@updateFilter);
    t2 = uicontrol(f, 'Style','text','Units','normalized',...
                   'Position',[0.91 0.32 0.07 0.04],...
                   'String', sprintf('%.0f Hz', default_high_cut));

    % --- Play Button ---
    uicontrol(f, 'Style','pushbutton','Units','normalized',...
              'String','Play Filtered Audio',...
              'Position',[0.35 0.22 0.3 0.05],...
              'Callback',@playFiltered);

    % --- Reset to Default Button ---
    uicontrol(f, 'Style','pushbutton','Units','normalized',...
              'String','Reset to Default',...
              'Position',[0.35 0.15 0.3 0.05],...
              'Callback',@resetDefaults);

    % --- Human Speech Range Info ---
    uicontrol(f, 'Style','text','Units','normalized',...
              'Position',[0.2 0.06 0.6 0.04],...
              'FontSize', 10,...
              'ForegroundColor', [0 0 0.5],...
              'String','Typical human speech frequency range: 300 Hz – 3400 Hz');

    % Initial output signal
    filtered_audio = x;

    %% --- Update Filter Based on Slider Values ---
    function updateFilter(~, ~)
        low_cut = s1.Value;
        high_cut = s2.Value;

        t1.String = sprintf('%.0f Hz', low_cut);
        t2.String = sprintf('%.0f Hz', high_cut);

        if high_cut <= low_cut + 50
        %%
        % 
        %   for x = 1:10
        %       disp(x)
        %   end
        % 
            high_cut = low_cut + 50;
            if high_cut >= fs/2 - 1
                high_cut = fs/2 - 1;
            end
            s2.Value = high_cut;
            t2.String = sprintf('%.0f Hz', high_cut);
        end

        Wn = [low_cut high_cut] / (fs/2);
        Wn = max(min(Wn, 0.999), 0.001);

        try
            [b, a] = butter(4, Wn, 'bandpass');
            filtered_audio = filter(b, a, x);
            set(h_plot, 'YData', filtered_audio);
            drawnow;
        catch err
            disp('Error designing filter:');
            disp(err.message);
        end
    end

    %% --- Play the Filtered Audio ---
    function playFiltered(~, ~)
        sound(filtered_audio, fs);
    end

    %% --- Reset Sliders to Default and Update Filter ---
    function resetDefaults(~, ~)
        s1.Value = default_low_cut;
        s2.Value = default_high_cut;
        t1.String = sprintf('%.0f Hz', default_low_cut);
        t2.String = sprintf('%.0f Hz', default_high_cut);
        updateFilter();
    end
end

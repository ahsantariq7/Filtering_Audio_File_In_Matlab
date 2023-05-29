[file, path] = uigetfile({'*.mp3';'*.wav';'*.flac';'*.m4a';'*.aac';'*.ogg';'*.wma'}, 'Select an audio file');
if isequal(file,0)
   disp('User selected Cancel');
else
   % Read the audio file
   [s, fs] = audioread(fullfile(path, file));

   % Display information about the audio file
   disp(['File selected: ' file]);
   disp(['Sampling frequency: ' num2str(fs) ' Hz']);
   disp(['Duration: ' num2str(size(s,1)/fs) ' seconds']);

end


% Define analysis parameters
N = 1000; % FFT size
M = 1000; % Window size
win = hamming(M); % Analysis window
L = length(s); % Length of speech signal
n = 1:M:L-M+1; % Starting indices of analysis windows

% Initialize output signal
s_mod = zeros(size(s));
w = hamming(N);
% Loop over all analysis windows
for i = 1:length(n)
    
    % Extract current analysis window
    s_win = s(n(i):n(i)+M-1);
    
    % Apply analysis window
    s_win = s_win .* win;
    
    % Compute magnitude spectrum of windowed signal
    P = abs(fft(s_win, N)).^2 / N;
    
    % Modify spectrum to enhance speech quality
    % TODO: Modify P to enhance speech quality by modifying its spectrum
    
    % Synthesize modified speech signal
    s_mod_win = real(ifft(sqrt(P) .* exp(1i * angle(fft(s_win, N)))));
    disp(size(s_mod_win));
    disp(size(win));
    % Apply synthesis window
    s_mod_win = s_mod_win .* win;

    % Add windowed, synthesized speech to output signal
    s_mod(n(i):n(i)+M-1) = s_mod(n(i):n(i)+M-1) + s_mod_win;
end
% Plot the input speech signal, analysis and synthesis windows, and output speech signal
figure;
subplot(2,2,1);
plot(s);
title('Input Speech Signal');
xlabel('Time (samples)');
ylabel('Amplitude');
xlim([0, length(s)]);

subplot(2,2,2);
plot(win);
title('Analysis Window');
xlabel('Sample Number');
ylabel('Window Value');
xlim([0, M]);

subplot(2,2,3);
plot(w);
title('Synthesis Window');
xlabel('Sample Number');
ylabel('Window Value');
xlim([0, N]);

subplot(2,2,4);
plot(s_mod_win);
title('Output Speech Signal');
xlabel('Time (samples)');
ylabel('Amplitude');
xlim([0, length(s_mod_win)]);

% Normalize output signal
s_mod = s_mod / max(abs(s_mod));

% Write output signal to file
audiowrite('H:\matlab_tasks\speech_signal_mod1.wav', s_mod, fs);

% Load the modified speech signal from the file
[s_mod, fs] = audioread('speech_signal_mod1.wav');

% Plot the modified speech signal
figure;
plot(s_mod);
title('Modified Speech Signal');
xlabel('Time (samples)');
ylabel('Amplitude');


% Compute the power of the original speech signal
s_pow = mean(s.^2);

% Compute the power of the difference between the original and modified speech signals
diff_pow = mean((s - s_mod).^2);

% Compute the SNR in dB
snr = 10*log10(s_pow/diff_pow);

% Display the SNR
fprintf('SNR: %.2f dB\n', snr);

% Plot the SNR difference
figure;
bar([s_pow, diff_pow]);
xticklabels({'Original Signal', 'Modified Signal'});
ylabel('SNR (dB)');
legend('Before', 'After');
title('SNR Comparison');
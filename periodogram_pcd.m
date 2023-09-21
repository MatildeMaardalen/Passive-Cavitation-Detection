close all

%% 1. Setup

% 1.1 Paths
dataset_dir = "E:\PCD\September 13 2023\OVA2\OVA_093Vpn\";          % folder with .mat files

files = dir(dataset_dir);
files = files(~ismember({files.name}, {'.', '..'}));                % remove aliases for parent and current folder
n_files = length(files);

% 1.2 Parameters
pre_gain = 5;                                                       % preamplifier gain (multiplication)
f0 = 0.5e6;                                                         % centre frequency of HIFU transducer (Hz)              
window = 2000:3500;                                                 % (optional) pulse window
background_win = [0,3];                                             % time window [s] of background signal
signal_win = [5,20];                                                % time window [s] of actual signal


%% 2. Loop over all files
% 2.1 Initialise variables
time = zeros(n_files,1);  
PSD = struct();
    
   
% 2.2 Loop over files
for i = 1:n_files
        
    % 2.2.1 Load file
    load(fullfile(files(i).folder,files(i).name));
    
    if i == 1
        fs = tpd.SampleFrequency; 
        [~,F] = pwelch(detrend(double(tpd.Data(window)))/pre_gain,[],[],[],fs);
        PSD.total       = zeros(length(F),n_files);                           % total PSD    
    end
   
        
    % 2.2.2 Get time point of file
    tmp = datevec(tpd.DateTime);
    time(i) = tmp(4)*60^2 + tmp(5)*60 + tmp(6);                 % hour, min, sec
        
    % 2.2.3 Power spectral density
    [PSD.total(:,i),~] = pwelch(double(tpd.Data(window))/pre_gain,[],[],[],fs); 
end

%% 3. Process PSD data
    
% 3.1 Sort by file acquisition time
[file_times,ind] = sort(time);
file_times       = file_times - file_times(1);                  % set time relative to first acquisiton
    
% 3.2 Sort the PSDs
PSD.total          = PSD.total(:, ind);

% 3.3 Extract background signal
background_win = file_times>background_win(1)...                % logic of file_times
   & file_times<background_win(2);

PSD.background = PSD.total(:, background_win);
PSD_background = mean(PSD.background,2);                        % calculate the average background signal

% 3.4 Extract actual signal
signal_win = file_times>signal_win(1)...                        % logic of file_times
   & file_times<signal_win(2);

PSD.signal = PSD.total(:, signal_win);
PSD_signal = mean(PSD.signal,2);                                % calculate the average signal

% 3.5 Find fmax limit (10 db above background)
diff = 0.5*db(PSD_signal)-0.5*db(PSD_background);               
[~,ind] = min(abs(diff-10));                                    % fmax should be where signal is ~10 db re V^2/Hz above background
fmax = F(ind)*1e-6;

%% 4. Plot signal and background PSDs

% 4.1 Create figure
figure;
plot(F*1e-6,0.5*db(PSD_background),'b-','LineWidth',1);
hold on;
plot(F*1e-6,0.5*db(PSD_signal),'r-','LineWidth',1);
grid on;
xlim([0 fs/2*1e-6])
fmin = 3;                                                       % found by looking at the plot
xline(fmin,'--b',{'fmin'})
xline(fmax,'--b',{'fmax'})
xlabel('Frequency [MHz]')
ylabel('PSD [dB re V^2/Hz]')




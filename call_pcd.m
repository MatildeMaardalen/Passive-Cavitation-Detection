
close all

%% 1. Setup

% 1.1 Paths
dataset_dir = "E:\PCD\2023\September 14 2023\MSA2\";

folders = dir(dataset_dir);
folders = folders(~ismember({folders.name}, {'.', '..'}));      % remove aliases for parent and current folder
folders = folders([folders(:).isdir]);                          % only list folders and not files in dataset_dir

% 1.2 Parameters
pre_gain = 5;                                                   % preamplifier gain (multiplication)
f0 = 0.5e6;                                                     % centre frequency of HIFU transducer (Hz)              
win = [2000,3500];                                              % (optional) pulse window
rec_win = [4,26];                                               % (optional) window to remove files with no US signal
fmin = 3e6;                                                     % exclude signal with harmonics (Hz)
fmax = 10.7e6;                                                  % exclude signal closer than 10db re V^2/Hz to background level

%% 2. Loop over all folders

% 2.1 Initialise variables
energy = zeros(length(folders),4);   
row_names = zeros(length(folders),1);
sample_info = {};

for i = 1:length(folders)
    
    % 2.2 Check if folder is empty
    folder_path = fullfile(dataset_dir,folders(i).name,'\','*.mat');
    if isempty(dir(folder_path))
        continue; 
    end
    
    % 2.3 Get parameteres
    folder_name = folders(i).name;
    sample_info = split(folder_name,'_');       
    voltage = str2double(sample_info{2}(1:3));                     % get voltage from file name
    pressure = 0.02472*voltage+0.09415;                            % pressure (PNP)-voltage (PNV) calibration
    
    file_name = strcat(sample_info{1}, " ",...
        num2str(round(pressure,1)), " MPa");
    
    % 2.4 Get power and energy data
    [Pwr, E] = fprocess_pcd(dataset_dir, folder_path, pre_gain,...
        f0, file_name, fmin, fmax, rec_win, win);
    
    % 2.5 Save energy data for each pressure in matrix
    energy(i,:) = E;
    row_names(i,1) = pressure;
    
end

%% Export energy data in csv file 
energy = [round(row_names,1) energy];
file_name = strcat('Energy_',sample_info{1});

results_table = array2table([energy],...
    'VariableNames', {'Pressures', 'Total Electrical Energy',...
    'Total Harmonic Electrical Energy',...
    'Total Ultraharmonic Electrical Energy',...
    'Total Broadband Electrical Energy'});
writetable(results_table, fullfile(dataset_dir,file_name));

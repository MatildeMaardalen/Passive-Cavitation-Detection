
function [Pwr, E] = fprocess_pcd(dataset_dir, folder_path, pre_gain,...
        f0, file_name, fmin, fmax, rec_win, win)
    
    %% 1. Setup
    
    % 1.1 Paths
    files = dir(folder_path);
    
    % 1.2 Load first file to get processing parameters
    load(fullfile(files(1).folder,files(1).name));
    fs = double(tpd.SampleFrequency);                               % sampling frequency
    
    % 1.3 Set processing window  
    if nargin==9
        if length(win)==1           
            win(2) = length(tpd.Data);
        end
    else        
        win(1) = double(tpd.PreSampleCount) + 1;                    
        win(2) = length(tpd.Data);                                
    end
   
    window = win(1):win(2); 
    
    % 1.4 Define variables
    % 1.4.1 Frequency parameters
    % fmin = 2*f0;                                                    % exclude fundamental                             
    % fmax = 0.9*(fs / 2);                                            % highest frequency to use (fs/2 is the nyquist frequency)
    fh =  fmin:f0:fmax;                                             % harmonics
    fuh = (fmin+0.5*f0):f0:fmax;                                    % ultraharmonics
    lfh = min([length(fh) length(fuh)]);
    
    % 1.4.2 Spectrum parameters
    bin = f0/20;                                                    % spectrum resolution (Hz)
    Nfft = round(fs/bin);                                           % samples needed for desired binsize
    Novr = 0.5*Nfft;                                                % sample overlap Welch method 
    
    % 1.4.3 Define bin width to get harmonics and ultraharmonics
    numbin_h = round(5.5e4/bin)-1;                                  % frequency range necessary to capture harmonic
    numbin_uh = round(2.5e4/bin);                                   % frequency range necessary to capture ultraharmonic
    numbin = [numbin_h numbin_uh];                                  % number of bins around each tone line [harm, ult. harm]
    bw{1} = -(numbin(1)-1)/2:(numbin(1)-1)/2;                       % relative bin indices for harmonic
    bw{2} = -(numbin(2)-1)/2:(numbin(2)-1)/2;                       % relative bin indices for ultraharmonic
    
    % 1.4.4 Initialise indice variable for harmonics and ultraharmonics
    ifh = zeros(numbin(1)*lfh,1);               
    ifu = zeros(numbin(2)*lfh,1);               

    % 1.5 Run first spectrum to confirm data size
    [~,F] = pwelch(detrend(double(tpd.Data(window)))/pre_gain,Nfft,Novr,Nfft,fs); 
    
    % 1.6 Get indices for harmonics and ultraharmonics
    for i = 1:lfh
       [~,ind] = min(abs(F-fh(i)));                                 % finds index where F is a harmonic
       ifh((1:numbin(1))+(i-1)*numbin(1)) = ind + bw{1};            % store the surrounding indices
       
       [~,ind] = min(abs(F-fuh(i)));                                % finds index where F is an ultraharmonic
       ifuh((1:numbin(2))+(i-1)*numbin(2)) = ind + bw{2};           % store the surrounding indices
    end
    
    % 1.7 Get frequency [MHz] values of indices
    f_MHz   = F*1e-6;                                               % convert to MHz
    fh_MHz  = f_MHz(ifh);
    fuh_MHz = f_MHz(ifuh);
    
    % 1.8 Initialise data storage 
    n_files = length(files);
    
    PSD.total = zeros(length(F),n_files);                           % total PSD    
    PSD.harm  = zeros(length(ifh), n_files);                        % harmonic PSD
    PSD.ultra = zeros(length(ifuh), n_files);                       % ultraharmonic PSD
    
    %% 2. Loop over all files
    % 2.1 Initialise time variable
    time = zeros(n_files,1);                                                                   
    
    % 2.2 Progressbar
    wb = waitbar(0,'file loop');
    
    % 2.3 Loop over files
    for i = 1:n_files
        
        % 2.3.1 Load file
        load(fullfile(files(i).folder,files(i).name));
        
        % 2.3.2 Get time point of file
        tmp = datevec(tpd.DateTime);
        time(i) = tmp(4)*60^2 + tmp(5)*60 + tmp(6);                 % hour, min, sec
        
        % 2.3.3 Power spectral density
        
        [PSD.total(:,i),~] = pwelch(double(tpd.Data(window))/pre_gain,Nfft,Novr,Nfft,fs); 
        
        % 2.3.4 Extract harmonic and ultraharmonic PSDs
        PSD.harm(:,i)  = PSD.total(ifh,i);
        PSD.ultra(:,i) = PSD.total(ifuh,i);
        
        % 2.3.5 Update progressbar
        waitbar(i/n_files,wb)
    end
    close(wb)
    
    %% 3. Process PSD data
    
    % 3.1 Sort by file acquisition time
    [file_times,ind] = sort(time);
    file_times       = file_times - file_times(1);                  % set time relative to first acquisiton
    
    % 3.2 Sort the PSDs
    PSD.total          = PSD.total(:, ind);
    PSD.harm           = PSD.harm(:, ind);
    PSD.ultra          = PSD.ultra(:, ind);
    
    % 3.3 Remove signal recorded before US was turned on
    rec_win     = file_times>rec_win(1)...                          % logic rec_win of file_times
        & file_times<rec_win(2);
    file_times  = file_times(rec_win);                         
    file_times  = file_times - file_times(1);                       % reset it relative to the first acquisition
    
    PSD.total   = PSD.total(:, rec_win);
    PSD.harm    = PSD.harm(:, rec_win);
    PSD.ultra   = PSD.ultra(:, rec_win);
    
    % 3.4 Convert PSD to power
    ze = 50;                                                        % impedance [Ohms]
%     ub = min([ifh(1) ifuh(1)]):max([ifh(end) ifuh(end)]);           % freq range for power calcs, spanning Ph, Pu        
    [~,ifmin] = min(abs(F-fmin));
    [~,ifmax] = min(abs(F-fmax));
    
    ub = ifmin:ifmax;                                               % freq range for power calcs
    Pwr(:,1) = sum(PSD.total(ub,:),1)*bin/ze;                       % total power: Pwr = V^2/R
    Pwr(:,2) = sum(PSD.harm,1)*bin/ze;                              % harmonic
    Pwr(:,3) = sum(PSD.ultra,1)*bin/ze;                             % ultraharmonic
    Pwr(:,4) = Pwr(:,1)-sum(Pwr(:,2:3),2);                          % broadband
    
    % 3.5 Calculate the energy
    time_win = diff(window([1,end]))/fs;                            % pulse window in seconds
    E = sum(Pwr,1)*time_win;                                        % Energy
    
    %% 4. Display
    
    mx = 0.5*db(max(PSD.total(:)));                                % define colourbar maximum
    
    % 4.1 Plot spectrogram
    figure
    set(gcf, 'name', 'Spectrogram','Position',[10 160 920 710])
    pcolor(file_times,f_MHz, 0.5*db(PSD.total))
    shading flat, colormap jet
    xlabel('Time [s]'), ylabel('Frequency [MHz]')                  % x and y axis labels
    ax = gca; ax.FontSize = 20;                                    %axis fontsize
    title(file_name)
    xlim(file_times([1 end])), ylim([0 fs*1e-6/2])                 % x and y axis limits
    hb = colorbar('vert'); hb.Label.String='PSD [dB reV^2/Hz]';    % colourbar
    set(hb, 'FontSize', 20); caxis([-180 -80])                     % colourbar limits
    hold on; 
    fmin = min([fh_MHz; fuh_MHz]);
    plot(file_times([1 end]), fmin*[1 1],'--w',...
        file_times([1 end]), fmax*[1 1]*1e-6,'--w')
    

    % 4.1.1 Save spectrogram
    fig_name = strcat(strrep(file_name,' ','_'),'_','spectrogram.png');
    saveas(gcf, fullfile(dataset_dir,fig_name))
    
    % 4.2 Plot spectral content
    figure
    set(gcf, 'name', 'Spectral Content','Position',[10 160 920 710])
    plot(file_times, 0.5*db(Pwr(:,2:4)))
    xlabel('Time [s]'), ylabel('Power [db re W]')                  % x and y axis labels
    legend({'Harmonic','Ultraharmonic','Broadband'},'location','best')
    ax = gca; ax.FontSize = 20;                                    %axis fontsize
    title(file_name)
    xlim(file_times([1 end]))                                      % x axis limits
         
    % 4.2.1 Save spectral content
    fig_name = strcat(strrep(file_name,' ','_'),'_','spectrogram.png');
    saveas(gcf, fullfile(dataset_dir,fig_name))
    
end


   
    
    
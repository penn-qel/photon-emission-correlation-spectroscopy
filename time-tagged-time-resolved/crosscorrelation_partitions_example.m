% Analysis of photon correlation data stored in TTTR format
% Note that these codes require routines in the tttr-functions directory

clear variables;

% NOTE: {} brackets are necessary
% Datapath can take multiple files

datapath = {'Data\1iPicoData_072121.mat'};

analysispath = 'tttr-functions'; % 'photon-emission-correlation-spectroscopy\time-tagged-time-resolved' < need functions stored here
matfilepath = 'Processed data/'; % to store results of analysis
addpath(analysispath)

%% Parse file names

rho = [];

for iFile = 1:length(datapath)
    dataFileStruct{iFile} = load(datapath{iFile});

    datadate = dataFileStruct{iFile}.name(end-9:end-4); % Data is stored in a directory named with this date

    idxTTTR = find(contains(dataFileStruct{iFile}.data,'tttr'));

    ptu_files{iFile} = dataFileStruct{iFile}.data(idxTTTR);

     %Check if background was recorded from tracking scan. If recorded,
     %use background and signal from tracking to calculate rho for
     %background subtraction
     BkgdIndex = find(contains(dataFileStruct{1}.sDataHeader,'background (kCts/s)'));
     if ~isempty(BkgdIndex)
         bkgd = dataFileStruct{iFile}.sData(:,BkgdIndex);
         signal = dataFileStruct{iFile}.sData(:,BkgdIndex-1);
         rho = [rho; 1 - bkgd./signal];
     else
         rho = [];
     end
  
end

% Calculate rho and rho error for background correction
rhoMean = mean(rho);
rhoErr = std(rho)/length(rho);

% Combine all cells of ptu_files in one
ptu_files = cat(2, ptu_files{:});

%% Setting up parameters for processing data
% USER NEEDS TO SET PARAMETERS HERE

tauOffset = 0.5; % (ns) offset to delay=0
tauLimitLin = 10; % (ns) limits on delay for linear processing
tauLimitLog = 1e6;% (ns) limits on delay for log processing
tauRes = 2; %(ns) delay bin width for linear processing
tRes = 0.01;  % (s) Resolution of count-rate calculation when partitioning data
threshold = 32; % intensity threshold for partitioning data

scaletype = 'linear'; %set scaletype to 'log' for log processing of data or 'linear' for linear processing of data


%% Set up parameters based on user entered parameters in previous section

parseOpts.integerTime = false;
parseOpts.verbose = true;
parseOpts.vectorized = true;

clear T2opts;

% To enable calculation of time-averaged count rates, set t_res:
T2opts.tRes = tRes; % (s) Resolution of count-rate calculation

% Status update flags
T2opts.verbose = true;
T2opts.statusbar = true;
nfiles = length(ptu_files);
T2opts.countRateRanges = 1e3*[[0,threshold];[threshold,Inf]]; % uncomment to partition data based on measured intensity

% To define a uniform vector of delay bins, set tau_res and tau_limits:
T2opts.tauRes = 1e-9*tauRes; % (s) delay bin width

if strcmp(scaletype,'log')
    % Define an arbitrary nonuniform vector of delay bins:
    delayVec = [1e-9*(0.1:0.1:0.9),logspace(-9,log10(tauLimitLog*1e-9),30)];%,logspace(-3.7,log10(tauLimitLog*1e-9),15)]; % (s) Vector of delay values to use symmetrically on both sides of zero
    T2opts.tauAxis = sort(1e-9*tauOffset + [-delayVec,delayVec]); % symmetric vector of non-uniformly spaced delay values
else
    T2opts.tauLimits = 1e-9*tauOffset + 1.0e-9*[-tauLimitLin tauLimitLin]; % (s) [Lower,Upper] limits on delay for linear processing
end
 
%% Process each individual file
tStart = tic;

for ff=1:nfiles
    
    TTTRdata = TTTR_import_PTU(ptu_files{ff});
    parsedData = TTTR_extract_channel_times(TTTRdata,parseOpts);
    
    fprintf(1,'Processing file %u of %u...',ff,nfiles);
    
    AllT2data(ff) = TTTR_cross_correlation(parsedData,T2opts);
    
    fprintf(1,'Done!\n');
end

toc(tStart);


%% Saving AllT2data and AllT2bkgd for future data processing


sampleName = 'example';

if strcmp(scaletype,'log')
    datafilename = [sampleName,'_AllT2data_logScale_',num2str(tauOffset),...
        'nsDelay_',num2str(1e3*tRes),'msTRes_','1e',num2str(log10(tauLimitLog)),'nsTauMax_',num2str(threshold),'kctsPart.mat'];
else
    datafilename = [sampleName,'_AllT2data_',num2str(1e3*tauRes),'psBin_',num2str(tauOffset),...
        'nsDelay_',num2str(1e3*tRes),'msTRes_',num2str(tauLimitLin),'nsTauMax_',num2str(threshold),'kctsPart.mat'];
end

save(fullfile(matfilepath,datafilename),'AllT2data');

%% Plot - Can either plot data just processed above, or can load previously-processed saved data

% load processed data if it's already saved
%linear data:
%load('Processed data\example_AllT2data_2000psBin_0.5nsDelay_10msTRes_10nsTauMax_32kctsPart.mat') 
%scaletype = 'linear';
%log data:
%load('Processed data\example_AllT2data_logScale_0.5nsDelay_10msTRes_1e6nsTauMax_32kctsPart.mat')
%scaletype = 'log';



% Merge individual sets
T2data = TTTR_merge_T2data(AllT2data);
tauOffset = 0.5;


% plot partitioned counts vs time
figure;
TotalCR = sum(T2data.countRates,1);
plot(T2data.tAxis(T2data.tFlag==1),TotalCR(T2data.tFlag==1)*1e-3,'.');%,'LineWidth',0.5,'MarkerEdgeColor','None','MarkerFaceColor','Black')
hold on;
plot(T2data.tAxis(T2data.tFlag==2),TotalCR(T2data.tFlag==2)*1e-3,'.');%,'LineWidth',0.5,'MarkerEdgeColor','None')
xlabel('t(s)');%,'Interpreter','Latex');
ylabel('I_{PL} (kcts/s)')
title('Example counts vs time data partitioned at 32kcts/s')


%plot autocorrelation

if strcmp(scaletype,'log')
    figure;
    % Determine delay bins to exclude due to detector afterflashes
    excludeRange1 = 1e-9*tauOffset+1e-9*[-44,-30];
    excludeRange2 = 1e-9*tauOffset+1e-9*[30,44];
    exclude = TTTR_exclude_bins(T2data.tauAxis,excludeRange1) |...
        TTTR_exclude_bins(T2data.tauAxis,excludeRange2);% |...

    % Fold symmetric results into 1-sided correlation function
    T2data = TTTR_fold_T2data(T2data,1e-9*tauOffset,exclude);

    errorbar(1e9*T2data.tauAxis,T2data.g2(:,1),T2data.g2Lerr(:,1),T2data.g2Uerr(:,1),'.')
    hold on;
    errorbar(1e9*T2data.tauAxis,T2data.g2(:,2),T2data.g2Lerr(:,2),T2data.g2Uerr(:,2),'.')
    set(gca,'XScale','log')
    xlabel('\tau (ns)');
    ylabel('g^{(2)}(\tau)');
    title('Example logscale autocorrelation data partitioned at 32kcts/s')

else
    figure; 
    errorbar(1e9*T2data.tauAxis-tauOffset,T2data.g2(:,1),T2data.g2Lerr(:,1),T2data.g2Uerr(:,1),'.')
    hold on;
    errorbar(1e9*T2data.tauAxis-tauOffset,T2data.g2(:,2),T2data.g2Lerr(:,2),T2data.g2Uerr(:,2),'.')
    xlabel('\tau (ns)');
    ylabel('g^{(2)}(\tau)');
    title('Example linear scale autocorrelation data partitioned at 32kcts/s')
end




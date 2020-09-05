function [t,PER_sum,PER_VEC] = Inside_Iter(Payload_bit, PL_dB,DS_ns, RX_Antenna,Coderate, BW_MHz, N_PRB,Subspacing_kHz,SymNum )

t0 = tic;
simParameters = [];             % Clear simParameters variable
simParameters.NFrames = 20;      % Number of 10ms frames

%%%%%%%%%% Constant%%%%%%%%%%%%%%%%%%
const_BOLTZMANN                 = 1.38064852e-23;
Temperature                     = 293;

%%%%%%%%%%%%%%%%Tunable%%%%%%%%%%%%%%%%%%%%%%%%%
payloadlength = Payload_bit;   %80; % PPDU length
Ptx = 20; % dBm
PL = PL_dB;%85; %dB PL_VEC = 85;
Iter_Num = 1e3;
Workers = 500;
PER_VEC = zeros(Iter_Num/Workers,1);
%SNR_VEC = zeros(Workers,1);




%%
% The variable |displaySimulationInformation| controls the display of
% simulation information such as the HARQ process ID used for each
% subframe. In case of CRC error, the value of the index to the RV sequence
% is also displayed.

displaySimulationInformation = false;

%% UE and PUSCH Configuration
% Set the key parameters of the simulation. These include:
% 
% * The bandwidth in resource blocks (12 subcarriers per resource block)
% * Subcarrier spacing: 15, 30, 60, 120, 240 (kHz)
% * Cyclic prefix length: normal or extended
% * Cell ID
% * Number of transmit and receive antennas
% 
% A substructure containing the UL-SCH and PUSCH parameters is also
% specified. This includes:
% 
% * Target code rate
% * Allocated resource blocks (PRBSet)
% * Modulation scheme: 'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'
% * Number of layers
% * Transform precoding (enable/disable)
% * PUSCH transmission scheme and MIMO precoding matrix indication (TPMI)
% * Number of antenna ports
% * PUSCH mapping type
% * DM-RS configuration parameters
% 
% Other simulation wide parameters are:
% 
% * Propagation channel model: 'TDL' or 'CDL'
%
% Note that if transform precoding is enabled, the number of layers should
% be set to 1.

% Bandwidth, numerology (SCS and CP type) and other general parameters
BW = BW_MHz * 1e6; %5e6;                              % Hz, ????????????
simParameters.NRB = N_PRB;   %11;                % Bandwidth in number of resource blocks (52RBs at 15kHz SCS for 10MHz BW)
simParameters.SubcarrierSpacing =  Subspacing_kHz; %30;  % 15, 30, 60, 120, 240 (kHz)
Occu_BW = simParameters.NRB*simParameters.SubcarrierSpacing*12*1e3;
%%%%%%%%%%%%%%%%%%
assert( 0.5*BW < Occu_BW &&  Occu_BW < BW,'The bandwidth and PRB number do not match')
%%%%%%%%%%%%%%%
simParameters.CyclicPrefix = 'Normal'; % 'Normal' or 'Extended'
simParameters.NCellID = 0;             % Cell identity
simParameters.NTxAnts = 1;             % Number of transmit antennas
simParameters.NRxAnts = RX_Antenna;             % Number of receive antennas

% UL-SCH/PUSCH parameters
simParameters.PUSCH.TargetCodeRate = Coderate;%490 / 1024;      % Code rate used to calculate transport block sizes
simParameters.PUSCH.PRBSet = (0:simParameters.NRB-1); % PUSCH PRB allocation
simParameters.PUSCH.SymbolSet = 0:SymNum-1;            % PUSCH symbol allocation in each slot 0:1
simParameters.PUSCH.NohPRB = 0;                  % Additional RE overhead per PRB
simParameters.PUSCH.EnableHARQ = false;           % Enable/disable HARQ, if disabled, single transmission with RV=0, i.e. no retransmissions
simParameters.PUSCH.Modulation = 'QPSK';         % 'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'
simParameters.PUSCH.NLayers = 1;                 % Number of PUSCH layers
simParameters.PUSCH.RNTI = 1;                    % Radio Network Temporary Identifier
simParameters.PUSCH.TransformPrecoding = false;  % Enable/disable transform precoding
simParameters.PUSCH.TxScheme = 'nonCodebook';    % Transmission scheme ('nonCodebook','codebook')
simParameters.PUSCH.NAntennaPorts = 1;           % Number of antenna ports for codebook based precoding
simParameters.PUSCH.TPMI = 0;                    % Precoding matrix indicator for codebook based precoding
% PUSCH DM-RS configuration
simParameters.PUSCH.DMRSSymbolSet = 0;
simParameters.PUSCH.PUSCHMappingType = 'B';      % PUSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
%simParameters.PUSCH.DMRSTypeAPosition = 2;       % Mapping type A only. First DM-RS symbol position (2,3)
simParameters.PUSCH.DMRSLength = 1;              % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
simParameters.PUSCH.DMRSAdditionalPosition = 0;  % Additional DM-RS symbol positions (max range 0...3)
simParameters.PUSCH.DMRSConfigurationType = 2;   % DM-RS configuration type (1,2)
%simParameters.PUSCH.NumCDMGroupsWithoutData = 2; % CDM groups without data
simParameters.PUSCH.NIDNSCID = 0;                % Scrambling identity (0...65535)
simParameters.PUSCH.NSCID = 0;                   % Scrambling initialization (0,1)
simParameters.PUSCH.NRSID = 0;                   % Scrambling ID for low-PAPR sequences (0...1007)
simParameters.PUSCH.GroupHopping = 'Disable';    % Hopping type ('Enable','Disable')

% Define the propagation channel type
simParameters.ChannelType = 'TDL'; % 'CDL' or 'TDL'

%%
% Create UE configuration structure |ue| and PUSCH configuration structure 
% |pusch|.

ue = simParameters;
pusch = simParameters.PUSCH;
pusch.NSlot = 0;

%%
% For key simulation parameters, define local variables for convenience.

%snrIn = simParameters.SNRIn;
nTxAnts = simParameters.NTxAnts;
nRxAnts = simParameters.NRxAnts;
channelType = simParameters.ChannelType;
%SymNum = length(simParameters.PUSCH.SymbolSet);

%% Propagation Channel Model Configuration
% Create the channel model object. Both CDL and TDL channel models are
% supported [ <#14 5> ].

if strcmpi(channelType,'CDL')
    
    channel = nrCDLChannel;
    channel.DelayProfile = 'CDL-A';
    [txsize,rxsize] = hArrayGeometry(nTxAnts,nRxAnts,'uplink');
    channel.TransmitAntennaArray.Size = txsize;
    channel.ReceiveAntennaArray.Size = rxsize;
    
else
    
    channel = nrTDLChannel;
    channel.DelayProfile = 'TDL-A';
    channel.NumTransmitAntennas = nTxAnts;
    channel.NumReceiveAntennas = nRxAnts;
    
end

channel.DelaySpread = DS_ns*1e-9;   %55e-9; % in seconds
freq = 3.5e9;
v = 2.5;
lambda = physconst('LightSpeed')/freq;
dopplershift = speed2dop(v,lambda);

channel.MaximumDopplerShift = 30; % in Hz

%%
% The sampling rate for the channel model is set using the value returned 
% from hOFDMInfo.

waveformInfo = hOFDMInfo(ue);
channel.SampleRate = waveformInfo.SamplingRate;

% calculate thermal noise power (k*T*B)
noisePower = const_BOLTZMANN * Temperature * channel.SampleRate;
P_RX =  Ptx - PL;   %dBm
P_RX = 10^((P_RX-30)/10);
SNR = P_RX/noisePower;
SNRdB = 10*log10(SNR);

%%
% Get the maximum number of delayed samples by a channel multipath
% component. This is calculated from the channel path with the largest
% delay and the implementation delay of the channel filter. This is
% required later to flush the channel filter to obtain the received signal.

chInfo = info(channel);
maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate));
maxChDelay = maxChDelay + chInfo.ChannelFilterDelay;

%% Processing Loop
% To determine the throughput at each SNR point, the PUSCH data is analyzed 
% per transmission instance using the following steps:
% 
% * _Update current HARQ process._ Check the CRC of the previous
% transmission for the given HARQ process. Determine whether a
% retransmission is required. If that is not the case generate new data.
% * _Generate resource grid._ Channel coding is performed by
% <docid:5g_ref#mw_sysobj_nrULSCH nrULSCH>. It operates on the input
% transport block provided. Internally, it keeps a copy of the transport
% block in case a retransmission is required. The coded bits are modulated
% by <docid:5g_ref#mw_function_nrPUSCH nrPUSCH>. Implementation-specific 
% MIMO precoding is applied to the resulting signal. Note that if 
% |TxScheme='codebook'|, codebook based MIMO precoding will already have
% been applied inside |nrPUSCH| and the implementation-specific MIMO 
% precoding is an additional stage of MIMO precoding. 
% * _Generate waveform._ The generated grid is then OFDM modulated.
% * _Model noisy channel._ The waveform is passed through a CDL or TDL
% fading channel. AWGN is added. The SNR for each layer is defined per RE
% and per receive antenna.
% * _Perform synchronization and OFDM demodulation._ Information returned
% by the channel is used for perfect synchronization. The synchronized
% signal is then OFDM demodulated.
% * _Perform perfect channel estimation._ Perfect channel estimation is
% used.
% * _Extract PUSCH and perform equalization._ The resource elements
% corresponding to the PUSCH allocation are extracted from the received
% OFDM resource grid and the channel estimate using
% <docid:5g_ref#mw_function_nrExtractResources nrExtractResources>. The
% received PUSCH resource elements are then MMSE equalized using
% <docid:5g_ref#mw_function_nrEqualizeMMSE nrEqualizeMMSE>.
% * _Decode the PUSCH._ The equalized PUSCH symbols, along with a noise
% estimate, are demodulated and descrambled by
% <docid:5g_ref#mw_function_nrPUSCHDecode nrPUSCHDecode> to obtain an
% estimate of the received codewords.
% * _Decode the Uplink Shared Channel (UL-SCH) and store the block CRC
% error for a HARQ process._ The vector of decoded soft bits is passed to
% <docid:5g_ref#mw_sysobj_nrULSCHDecoder nrULSCHDecoder> which decodes
% the codeword and returns the block CRC error used to determine the
% throughput of the system.



% Set up Redundancy Version (RV) sequence, number of HARQ processes and
% the sequence in which the HARQ processes are used
if pusch.EnableHARQ
    % From PUSCH demodulation requirements in RAN WG4 meeting #88bis
    % (R4-1814062)
    rvSeq = [0 2 3 1];
else
    % HARQ disabled - single transmission with RV=0, no retransmissions
    rvSeq = 0;
end
% Specify the order in which we cycle through the HARQ processes
% NHARQProcesses = 16;
% harqSequence = 1:NHARQProcesses;

 % Calculate the transport block size 
[puschIndices,dmrsIndices,dmrsSymbols,puschIndicesInfo] = hPUSCHResources(ue,pusch);
TBS = hPUSCHTBS(pusch,puschIndicesInfo.NREPerPRB - pusch.NohPRB);
assert(TBS>=payloadlength,'The symbol number is not enough to carry the PPDU bits')
Bit_Capacity = puschIndicesInfo.G;

% Create UL-SCH encoder System object
encodeULSCH = nrULSCH;
encodeULSCH.MultipleHARQProcesses = false;
encodeULSCH.TargetCodeRate = pusch.TargetCodeRate;

% Create UL-SCH decoder System object
decodeULSCH = nrULSCHDecoder;
decodeULSCH.MultipleHARQProcesses = false;
decodeULSCH.TargetCodeRate = pusch.TargetCodeRate;
decodeULSCH.TransportBlockLength = TBS;

% % Total number of OFDM symbols in the simulation period
% NSymbols = ue.NFrames * 10 * waveformInfo.SymbolsPerSubframe;
%     
% % OFDM symbol number associated with start of each PUSCH transmission
% ue.NSymbol = 0;




for iterIdx = 1:Iter_Num/Workers
    
    % Reset the random number generator and channel, so that each SNR point
    % will experience the same noise and channel realizations
%     rng('default');
%     reset(channel);
    
    % Initialize the state of all HARQ processes and reset the UL-SCH 
    % decoder
%     harqProcesses = hNewHARQProcesses(NHARQProcesses,rvSeq,1);
%     harqProcCntr = 0; % HARQ process counter
%     reset(decodeULSCH);
    
    %SNRdB = snrIn(snrIdx);
    %fprintf('\nSimulating %s-based transmission scheme (%dx%d) and SCS=%dkHz with %s channel at %gdB pathloss for %d 10ms frame(s)\n',pusch.TxScheme,nTxAnts,nRxAnts,ue.SubcarrierSpacing,channelType,PL,ue.NFrames);

    % Total number of OFDM symbols in the simulation period
%     NSymbols = ue.NFrames * 10 * waveformInfo.SymbolsPerSubframe;
%     
%     % OFDM symbol number associated with start of each PUSCH transmission
%     ue.NSymbol = 0;
    
    % Running counter of the number of PUSCH transmission instances
    % The simulation will use this counter as the slot number for each
    % PUSCH
    %pusch.NSlot = 0;
    %PER = 0;
    PER_worker_Vec = zeros(1,Workers);
    
    %while ue.NSymbol < NSymbols
    parfor PKT_Idx = 1:Workers
    
        trBlk = randi([0 1],TBS,1);
        setTransportBlock(encodeULSCH,trBlk);
     
        % UL-SCH encoding
        %codedTrBlock = encodeULSCH(pusch.Modulation,pusch.NLayers,puschIndicesInfo.G,harqProcesses(harqProcIdx).RV,harqProcIdx-1);
        %codedTrBlock = encodeULSCH(pusch.Modulation,pusch.NLayers,puschIndicesInfo.G,rvSeq);
        codedTrBlock = encodeULSCH(pusch.Modulation,pusch.NLayers,Bit_Capacity,rvSeq);


        % PUSCH modulation, including codebook based MIMO precoding if 
        % TxScheme = 'codebook'
        MRB = numel(pusch.PRBSet);
        puschSymbols = nrPUSCH(codedTrBlock,pusch.Modulation,pusch.NLayers,ue.NCellID,pusch.RNTI,pusch.TransformPrecoding,MRB,pusch.TxScheme,pusch.NAntennaPorts,pusch.TPMI);
        
        % Create resource grid associated with PUSCH transmission period
        puschGrid = zeros(waveformInfo.NSubcarriers,SymNum,nTxAnts);  %waveformInfo.SymbolsPerSlot
        %puschGrid = zeros(waveformInfo.NSubcarriers,SymNum,nTxAnts);
        % Implementation-specific PUSCH MIMO precoding and mapping. This 
        % MIMO precoding step is in addition to any codebook based 
        % MIMO precoding done during PUSCH modulation above
        if (strcmpi(pusch.TxScheme,'codebook'))
            % codebook based MIMO precoding, F precodes between PUSCH
            % transmit antenna ports and transmit antennas
            F = eye(pusch.NAntennaPorts,nTxAnts);
        else
            % non-codebook based MIMO precoding, F precodes between PUSCH 
            % layers and transmit antennas
            F = eye(pusch.NLayers,nTxAnts);
        end
        [~,puschAntIndices] = nrExtractResources(puschIndices,puschGrid);
        puschGrid(puschAntIndices) = puschSymbols * F;
        
        % Implementation-specific PUSCH DM-RS MIMO precoding and mapping.
        % The DM-RS creation in hPUSCHResources above includes codebook
        % based MIMO precoding if applicable
        for p = 1:size(dmrsSymbols,2)
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),puschGrid);
            puschGrid(dmrsAntIndices) = puschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * F(p,:);
        end
        
        % OFDM modulation
        txWaveform = hOFDMModulate(ue,puschGrid);
        txWaveform = txWaveform/norm(txWaveform);

        % Pass data through channel model. Append zeros at the end of the
        % transmitted waveform to flush channel content. These zeros take
        % into account any delay introduced in the channel. This is a mix
        % of multipath delay and implementation delay. This value may 
        % change depending on the sampling rate, delay profile and delay
        % spread
        txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))]; %#ok<AGROW>
        [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
        % Add AWGN to the received time domain waveform 
        % Normalize noise power by the IFFT size used in OFDM modulation,
        % as the OFDM modulator applies this normalization to the
        % transmitted waveform. Also normalize by the number of receive
        % antennas, as the default behaviour of the channel model is to
        % apply this normalization to the received waveform
        %SNR = 10^(SNRdB/20);
        
        
        P_rxWaveform = sum(var(rxWaveform)*size(rxWaveform,1));
        N0 = sqrt(P_rxWaveform/(2.0*SymNum*double(waveformInfo.Nfft)*SNR)); %    ,  *nRxAnts
        noise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
        rxWaveform = rxWaveform + noise;

        % Perfect synchronization. Use information provided by the channel
        % to find the strongest multipath component
        pathFilters = getPathFilters(channel);
        [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
        rxWaveform = rxWaveform(1+offset:end,:);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid
        rxGrid = hOFDMDemodulate(ue,rxWaveform);
        
        % Perfect channel estimation, use the value of the path gains
        % provided by the channel
        estChannelGrid = nrPerfectChannelEstimate_modi(pathGains,pathFilters,ue.NRB,ue.SubcarrierSpacing,pusch.NSlot,offset,sampleTimes,ue.CyclicPrefix);

        % Get perfect noise estimate (from the noise realization)
        noiseGrid = hOFDMDemodulate(ue,noise(1+offset:end,:));
        noiseEst = var(noiseGrid(:));
        
        % Apply MIMO precoding to estChannelGrid
        % Linearize 4D matrix and reshape after multiplication
        K = size(estChannelGrid,1);
        estChannelGrid = reshape(estChannelGrid,K*SymNum*nRxAnts,nTxAnts);  %waveformInfo.SymbolsPerSlot
        estChannelGrid = estChannelGrid * F.';
        if (strcmpi(pusch.TxScheme,'codebook'))
            W = nrPUSCHCodebook(pusch.NLayers,pusch.NAntennaPorts,pusch.TPMI,pusch.TransformPrecoding);
            estChannelGrid = estChannelGrid * W.';
        end
        estChannelGrid = reshape(estChannelGrid,K,SymNum,nRxAnts,[]); %waveformInfo.SymbolsPerSlot
        
        % Get PUSCH resource elements from the received grid
        [puschRx,puschHest] = nrExtractResources(puschIndices,rxGrid,estChannelGrid);
        
        % Equalization
        [puschEq,csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEst);
        
        % Decode PUSCH physical channel
        [ulschLLRs,rxSymbols] = nrPUSCHDecode(puschEq,pusch.Modulation,ue.NCellID,pusch.RNTI,noiseEst,pusch.TransformPrecoding,MRB);
        
        % Apply channel state information (CSI) produced by the equalizer,
        % including the effect of transform precoding if enabled
        if (pusch.TransformPrecoding)
            MSC = MRB * 12;
            csi = nrTransformDeprecode(csi,MRB) / sqrt(MSC);
            csi = repmat(csi((1:MSC:end).'),1,MSC).';
            csi = reshape(csi,size(rxSymbols));
        end
        csi = nrLayerDemap(csi);
        Qm = length(ulschLLRs) / length(rxSymbols);
        csi = reshape(repmat(csi{1}.',Qm,1),[],1);
        ulschLLRs = ulschLLRs .* csi;
        
        % Decode the UL-SCH transport channel
        %decodeULSCH.TransportBlockLength = TBS;
        [decbits] = decodeULSCH(ulschLLRs,pusch.Modulation,pusch.NLayers,rvSeq);
        PER_worker_Vec(PKT_Idx) = (~isequal(decbits(1:payloadlength),trBlk(1:payloadlength)));
        %PKT_VEC(PLIdx) =  PKT_VEC(PLIdx) + 1;
        % Store values to calculate throughput
%         simThroughput(snrIdx) = simThroughput(snrIdx) + (~harqProcesses(harqProcIdx).blkerr * TBS);
%         maxThroughput(snrIdx) = maxThroughput(snrIdx) + TBS;
        
        % Display transport block error information
%         if (displaySimulationInformation)
%             fprintf('\n(%3.2f%%) HARQ Proc %d: ',100*(ue.NSymbol+size(puschGrid,2))/NSymbols,harqProcIdx);
%             estrings = {'passed','failed'};
%             rvi = harqProcesses(harqProcIdx).RVIdx;
%             if rvi == 1
%                 ts = sprintf('Initial transmission (RV=%d)',rvSeq(rvi));
%             else
%                 ts = sprintf('Retransmission #%d (RV=%d)',rvi-1,rvSeq(rvi));
%             end
%             fprintf('%s %s. ',ts,estrings{1+harqProcesses(harqProcIdx).blkerr});
%         end
        
%         % Update starting symbol number of next PUSCH transmission
%         ue.NSymbol = ue.NSymbol + waveformInfo.SymbolsPerSlot;
%         
%         % Update count of overall number of PUSCH transmissions
%         pusch.NSlot = pusch.NSlot + 1;
        
        % Update HARQ process counter
        %harqProcCntr = harqProcCntr + 1;

    end
    PER_VEC(iterIdx) = sum(PER_worker_Vec);
    

end

%% Results
% Display the measured throughput. This is calculated as the percentage of
% the maximum possible throughput of the link given the available resources
% for data transmission.
%PER_VEC = PER_VEC/(ue.NFrames * 10* waveformInfo.SlotsPerSubframe);
PER_sum = sum(PER_VEC)/Iter_Num;
% Display the results dynamically in the command window
% if (displaySimulationInformation)
%     fprintf('\n');
% end
fprintf([['\nPER for ' num2str(Iter_Num) ' packet(s) '],'= %.10f\n'], PER_sum);%PER_VEC(PLIdx)/(ue.NFrames * 10* waveformInfo.SlotsPerSubframe))



t= toc(t0);
end





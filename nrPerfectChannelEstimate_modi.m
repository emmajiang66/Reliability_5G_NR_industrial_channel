function H = nrPerfectChannelEstimate_modi(pathGains,pathFilters,NRB,SCS,initialNSlot,varargin)
%nrPerfectChannelEstimate perfect channel estimation
%   H = nrPerfectChannelEstimate(...) performs perfect channel estimation,
%   producing a perfect channel estimate H, by reconstructing the channel
%   impulse response from information about the propagation channel and
%   then performing OFDM demodulation. H is a K-by-L-by-Nr-by-Nt array
%   where K is the number of subcarriers, L is the number of OFDM symbols,
%   Nr is the number of receive antennas and Nt is the number of transmit
%   antennas.
%
%   H = nrPerfectChannelEstimate(PATHGAINS,PATHFILTERS,NRB,SCS,INITIALNSLOT)
%   performs perfect channel estimation by reconstructing the channel
%   impulse response from the channel path gains array PATHGAINS and path
%   filter impulse response matrix PATHFILTERS, and then performing OFDM
%   demodulation for NRB resource blocks with subcarrier spacing SCS and
%   initial slot number INITIALNSLOT.
%
%   PATHGAINS must be an array of size Ncs-by-Np-by-Nt-by-Nr, where Ncs is
%   the number of channel snapshots and Np is the number of paths. The
%   times of the channel snapshots are given by the SAMPLETIMES input (see
%   below).
%
%   PATHFILTERS must be a matrix of size Nh-by-Np where Nh is the number of
%   impulse response samples.
%
%   NRB is the number of resource blocks (1...275).
%
%   SCS is the subcarrier spacing in kHz (15, 30, 60, 120, 240).
%
%   INITIALNSLOT is the 0-based initial slot number, a non-negative scalar
%   integer. The initial slot number is used modulo the number of slots per
%   subframe to select the appropriate cyclic prefix lengths for OFDM
%   demodulation.
%
%   H = nrPerfectChannelEstimate(...,OFFSET) allows the timing offset
%   OFFSET to be specified, an integer number of samples indicating where
%   the OFDM demodulation will start on the reconstructed waveform. The
%   default is to use <a href="matlab:doc('nrPerfectTimingEstimate')">nrPerfectTimingEstimate</a> to establish the timing
%   offset.
%
%   H = nrPerfectChannelEstimate(...,OFFSET,SAMPLETIMES) allows the sample
%   times SAMPLETIMES of the channel snapshots to be specified. SAMPLETIMES
%   must be of size Ncs-by-1 and specifies the time of occurrence of each
%   channel snapshot (the 1st dimension of PATHGAINS). The default is a
%   vector of times starting at zero, where the number of snapshots is
%   given by the 1st dimension sizing of PATHGAINS and the sampling rate is
%   equal to the sampling rate used for OFDM modulation for number of
%   resource blocks NRB and subcarrier spacing SCS.
%
%   H = nrPerfectChannelEstimate(...,CPLEN) allows the cyclic prefix length
%   to be specified. CPLEN must be 'normal' for normal cyclic prefix length
%   (default) or 'extended' for extended cyclic prefix length. Note that
%   for the numerologies specified in TS 38.211 Section 4.2, extended
%   cyclic prefix length is only applicable for 60 kHz subcarrier spacing.
%
%   % Example:
%   % Configure a TDL-C channel with 100ns delay spread and plot the 
%   % estimated channel magnitude response for the first receive antenna.
%
%   NRB = 25;
%   SCS = 15;
%   nSlot = 0;
%   SR = 7.68e6;
%   
%   tdl = nrTDLChannel;
%   tdl.DelayProfile = 'TDL-C';
%   tdl.DelaySpread = 100e-9;
%   tdl.MaximumDopplerShift = 300;
%   tdl.SampleRate = SR;
%   
%   T = SR * 1e-3;
%   tdlInfo = info(tdl);
%   Nt = tdlInfo.NumTransmitAntennas;
%   in = complex(randn(T,Nt),randn(T,Nt));
%
%   [~,pathGains] = tdl(in);
%   pathFilters = getPathFilters(tdl);
%
%   hest = nrPerfectChannelEstimate(pathGains,pathFilters,NRB,SCS,nSlot);
%   size(hest)
%
%   figure;
%   surf(abs(hest(:,:,1)));
%   shading('flat');
%   xlabel('OFDM symbols');
%   ylabel('Subcarriers');
%   zlabel('|H|');
%   title('Channel magnitude response');
%
%   % 2:
%   % Repeat the channel estimate for extended cyclic prefix.
%
%   hest = nrPerfectChannelEstimate(pathGains,pathFilters,NRB,SCS, ...
%           nSlot,'extended');
%   size(hest)
%
%   % 3:
%   % Configure a CDL-D channel with 30ns delay spread and plot the 
%   % estimated channel magnitude response for the first receive antenna.
%   
%   cdl = nrCDLChannel;
%   cdl.DelayProfile = 'CDL-D';
%   cdl.DelaySpread = 30e-9;
%   cdl.MaximumDopplerShift = 5;
%
%   cdlInfo = info(cdl);
%   Nt = cdlInfo.NumTransmitAntennas;
%   in = complex(randn(T,Nt),randn(T,Nt));
%
%   [~,pathGains,sampleTimes] = cdl(in);
%   pathFilters = getPathFilters(cdl);
%
%   offset = nrPerfectTimingEstimate(pathGains,pathFilters);
%
%   hest = nrPerfectChannelEstimate(pathGains,pathFilters,...
%              NRB,SCS,nSlot,offset,sampleTimes);
%
%   figure;
%   surf(abs(hest(:,:,1)));
%   shading('flat');
%   xlabel('OFDM symbols');
%   ylabel('Subcarriers');
%   zlabel('|H|');
%   title('Channel magnitude response');
%
%   See also nrPerfectTimingEstimate, nrTDLChannel, nrCDLChannel.

%   Copyright 2018 The MathWorks, Inc.

%#codegen

    narginchk(5,8);
    
    % Validate mandatory inputs
    validateInputs(pathGains,pathFilters,NRB,SCS,initialNSlot);
    
    % Parse and validate optional inputs, and calculate OFDM information 
    % structure 'ofdminfo'
    [offset,sampleTimes,ofdminfo] = getOptionalInputs(pathGains, ...
        pathFilters,NRB,SCS,varargin{:});
    
    % Get number of channel impulse response samples 'Nh'
    Nh = size(pathFilters,1);
    
    % Get number of paths 'Np', number of transmit antennas 'Nt' and number
    % of receive antennas 'Nr' in the path gains array. The pathGains are
    % of size Ncs-by-Np-by-Nt-by-Nr, where 'Ncs' is the number of channel
    % snapshots
    [~,Np,Nt,Nr] = size(pathGains);
    
    % Set the origin of the sample times to zero, and establish the total
    % duration 'T' in samples
    sampleTimes = sampleTimes - sampleTimes(1);
    T = round(sampleTimes(end) * ofdminfo.SampleRate) + 1;
    
    % Establish the total number of subframes spanned by 'T', rounded up,
    % which determines the required number of repetitions of the cyclic
    % prefix lengths
    samplesPerSubframe = ofdminfo.SampleRate * 1e-3;
    nSubframes = ceil(T / samplesPerSubframe);
    
    % Establish the starting and ending sample indices of each OFDM symbol
    % across the total number of subframes, taking into consideration the
    % initial slot number, and update the cyclic prefix lengths to span all
    % subframes
    cpLengths = repmat(ofdminfo.CyclicPrefixLengths,1,nSubframes);
    initialNSlot = mod(initialNSlot,ofdminfo.SlotsPerSubframe);
    cpLengths = circshift(cpLengths,-initialNSlot*ofdminfo.SymbolsPerSlot);
    fftLengths = [0 repmat(ofdminfo.Nfft,1,numel(cpLengths)-1)];
    symbolStarts = cumsum(cpLengths + fftLengths);
    symbolEnds = symbolStarts + ofdminfo.Nfft;
    ofdminfo.CyclicPrefixLengths = cpLengths;
    
    % Ensure that total duration 'T' is at least one slot
    T = min(T,symbolEnds(ofdminfo.SymbolsPerSlot));
    
    % Establish how many OFDM symbols 'L' are spanned by 'T' time samples
    % and round down to nearest whole slot
    L = numel(find(symbolEnds <= T));
%     L = L - mod(L,ofdminfo.SymbolsPerSlot);
    symbolStarts = symbolStarts(1:L);
    
    % Establish which OFDM symbol start times correspond to which channel
    % coefficient sample times. 'idx' is a vector of length 'L' indicating
    % the 1st dimension index of 'pathGains' for each OFDM symbol start
    % time
    if (~all(diff(sampleTimes)==ofdminfo.SampleRate))
        t = [0; sampleTimes + mean(diff(sampleTimes))/2];
    else
        t = sampleTimes;
    end
    symbolStartTimes = (symbolStarts + offset) / ofdminfo.SampleRate;
    idx = sum(bsxfun(@(x,y)x>=y,symbolStartTimes,t),1);
    
    % Prepare the path gains matrix by indexing using 'idx' to select a
    % first dimension element for each OFDM symbol start, and permute to
    % put the multipath components in the first dimension and switch the
    % antenna dimensions. The pathGains are now of size Np-by-Nr-by-Nt-by-L
    pathGains = pathGains(idx,:,:,:);
    pathGains = permute(pathGains,[2 4 3 1]);

    % Create channel impulse response array 'h' for each impulse response
    % sample, receive antenna, transmit antenna and OFDM symbol
    h = zeros(Nh,Nr,Nt,L,'like',pathGains);

    % For each path, add its contribution to the channel impulse response
    % across all transmit antennas, receive antennas and OFDM symbols
    for np = 1:Np

        h = h + bsxfun(@times,pathFilters(:,np),pathGains(np,:,:,:));

    end

    % Create the empty received waveform (for each transmit antenna)
    rxWave = zeros([T Nr Nt],'like',pathGains);
    
    % For each OFDM symbol, add the corresponding impulse response samples
    % across all transmit antennas and receive antennas to the received
    % waveform. Note that the impulse responses are positioned according to
    % the timing offset 'offset' and the channel filter delay so that
    % channel estimate produced is as similar as possible to that produced
    % for a filtered waveform (without incurring the time cost of the full
    % filtering)
    for l = 1:L

        tl = symbolStarts(l) - offset + (1:Nh);
        rxWave(tl,:,:) = rxWave(tl,:,:) + h(:,:,:,l);

    end

    % Remove any samples from the end of the received waveforms that
    % correspond to incomplete OFDM symbols
    rxWave = rxWave(1:symbolEnds(L),:,:);

    % For each transmit antenna, OFDM demodulate the received waveform
    % across all receive antennas to form the overall channel estimate
    % array
    H = zeros(ofdminfo.K,L,Nr,Nt,'like',pathGains);
    for nt = 1:Nt

        rxGrid = ofdmDemodulate(rxWave(:,:,nt),ofdminfo,L);
        H(:,:,:,nt) = rxGrid(:,:,1:Nr);

    end

end

function out = ofdmInfo(NRB,SCS,ECP)
    
    % Subcarrier spacing configuration
    % TS 38.211 Section 4.3.2
    mu = log2(SCS / 15);
    
    % Total number of subcarriers in the resource grid
    % TS 38.211 Section 4.4.4.1
    K = NRB * 12;
    
    % Choose IDFT size based on three rules:
    %   * power of 2 size
    %   * maximum occupancy (K / Nfft) of 85%
    %   * minimum of 128 (so that cyclic prefix lengths are always integer)
    nfft = max(power(2,ceil(log2(K / 0.85))),128);
    
    % OFDM sample rate given by reference subcarrier spacing (15kHz),
    % subcarrier spacing configuration 'mu' and IDFT size 'Nfft'
    SR = 2^mu * 15e3 * nfft;
    
    % OFDM symbols per slot and per subframe
    % TS 38.211 Section 4.3.2
    if ECP
        symbolsPerSlot = 12;
    else
        symbolsPerSlot = 14;
    end
    L = symbolsPerSlot * 2^mu;    
    
    % Nominal cyclic prefix lengths in LTE numerology (SCS=15kHz,Nfft=2048)
    % TS 38.211 Section 5.3.1
    if ECP
        NCP = 512;
    else
        NCP = 144;
    end
    
    % Adjust for IDFT size in use
    NCP = NCP / 2048 * nfft;
    
    % Create vector of cyclic prefix lengths across a subframe and adjust 
    % cyclic prefix lengths at start of each half subframe
    cpLengths = NCP * ones(1,L);
    cpLengths([1 1+L/2]) = (SR*1e-3 - (L*nfft + (L-2)*NCP)) / 2;    
    
    out.K = K;
    out.Nfft = nfft;
    out.SampleRate = SR;
    out.SymbolsPerSlot = symbolsPerSlot;
    out.SymbolsPerSubframe = L;
    coder.varsize('out.CyclicPrefixLengths',[1 Inf],[0 1]);
    out.CyclicPrefixLengths = cpLengths;
    out.SlotsPerSubframe = 2^mu;
    
end

function grid = ofdmDemodulate(waveform,ofdminfo,L)
    
    cpFraction = 0.55;
    symOffset = fix(ofdminfo.CyclicPrefixLengths(1:L) * cpFraction);
    
    firstSC = (ofdminfo.Nfft/2) - (ofdminfo.K/2) + 1;
    nullIndices = [1:(firstSC-1) (firstSC+ofdminfo.K):ofdminfo.Nfft].';
    
    grid = ofdmdemod(waveform,ofdminfo.Nfft, ... 
        ofdminfo.CyclicPrefixLengths(1:L),symOffset,nullIndices);

end

function validateInputs(pathGains,pathFilters,NRB,SCS,initialNSlot)
% Check inputs

    fcnName = 'nrPerfectChannelEstimate';
    
    % Validate channel path gains
    validateattributes(pathGains,{'double','single'}, ...
        {},fcnName,'PATHGAINS');
    coder.internal.errorIf(ndims(pathGains)>4, ...
        'nr5g:nrPerfectChannelEstimate:InvalidPathDims',ndims(pathGains));
    
    % Validate path filters impulse response
    validateattributes(pathFilters,{'double'}, ...
        {'real','2d'},fcnName,'PATHFILTERS');
    coder.internal.errorIf(size(pathGains,2)~=size(pathFilters,2), ...
        'nr5g:nrPerfectChannelEstimate:InconsistentPaths', ...
        size(pathGains,2),size(pathFilters,2));
    
    % Validate the number of resource blocks (1...275)
    validateattributes(NRB,{'numeric'}, ...
        {'real','integer','scalar','>=',1,'<=',275},fcnName,'NRB');
    
    % Validate subcarrier spacing input in kHz (15/30/60/120/240)
    validateattributes(SCS,{'numeric'}, ...
        {'real','integer','scalar'},fcnName,'SCS');
    validSCS = [15 30 60 120 240];
    coder.internal.errorIf(~any(SCS==validSCS), ...
        'nr5g:nrPerfectChannelEstimate:InvalidSCS', ...
        SCS,num2str(validSCS));
    
    % Validate zero-based initial slot number
    validateattributes(initialNSlot,{'numeric'}, ...
        {'real','nonnegative','scalar','integer'},fcnName,'INITIALNSLOT');
    
end

function [offset,sampleTimes,ofdminfo] = getOptionalInputs(pathGains,pathFilters,NRB,SCS,varargin)
    
    coder.extrinsic('nr5g.internal.parseOptions');
    narg = numel(varargin);
    fnname = 'nrPerfectChannelEstimate';
    
    % Find positions of optional arguments: offset, sampleTimes, and
    % char/string value for cyclic prefix length which can appear in any
    % position
    argpos = [0 0]; % positions of offset and sampleTimes, 0 if absent
    pos = 1;
    cparg = 0; % position of name of cyclic prefix length, 0 if absent    
    for i = 1:narg
        if (ischar(varargin{i}) || isstring(varargin{i}) || pos>2)
            if (cparg==0)
                cparg = i;
            end
        else
            argpos(pos) = i;
            pos = pos + 1;
        end
    end
    
    % Parse cyclic prefix length value or provide a default value
    if cparg~=0
        cplen = varargin{cparg};
        validateattributes(cplen,{'char' 'string'},{},fnname,'CPLEN');
        validatestring(cplen,{'normal','extended'},fnname,'CPLEN');
        ECP = startsWith('extended',cplen,'IgnoreCase',true);
    else
        ECP = false;
    end
    
    % Calculate OFDM information
    ofdminfo = ofdmInfo(NRB,SCS,ECP);
    
    % Validate offset or provide a default value
    if argpos(1)~=0
        offset = varargin{argpos(1)};
        Nh = size(pathFilters,1);
        validateattributes(offset,{'numeric'}, ...
            {'nonnegative','scalar'},fnname,'offset');
        coder.internal.errorIf(offset>(Nh-1), ...
            'nr5g:nrPerfectChannelEstimate:InvalidOffset', ...
            offset,Nh);
    else
        % Default: use nrPerfectTimingEstimate to establish the timing
        % offset
        offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    end
    
    % Validate sampleTimes or provide a default value
    if argpos(2)~=0
        sampleTimes = varargin{argpos(2)};
        Ncs = size(pathGains,1);
        validateattributes(sampleTimes,{'double'}, ...
            {'column','increasing'},fnname,'sampleTimes');
        coder.internal.errorIf(length(sampleTimes)~=Ncs, ...
            'nr5g:nrPerfectChannelEstimate:InvalidSampleTimes', ...
            length(sampleTimes),Ncs);
    else
        % Default: vector of times at the OFDM sampling rate, one for each
        % channel snapshot in 'pathGains' and starting at zero
        sampleTimes = (0:(size(pathGains,1)-1)).' / ofdminfo.SampleRate;
    end
    
end

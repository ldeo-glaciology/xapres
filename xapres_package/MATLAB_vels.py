import numpy as np
def MATLAB_vels.py(f,g):
    t1 = f.time.data
    t2 = g.time.data
    dt = (t2-t1)/ np.timedelta64(1, 's')
    dr = np.mean(np.diff(f.profile_range))
    # configurations
    maxStrain = 0.005
    maxDepth = 1400
    minDepth = 10
    coarseChunkWidth = 15
    
    # COARSE ALIGN
    maxOffset = maxStrain*(maxDepth-minDepth)
    stepSizeM = 5
    binStart = [minDepth:stepSizeM:maxDepth-coarseChunkWidth]
    AC_range = np.zeros(binStart.shape)
    AC_dh = np.zeros(binStart.shape)
    for ii in range(0,len(binStart)):
        depthRange = [binStart(ii) binStart(ii)+coarseChunkWidth];
        fi = np.argwhere((f.profile_range>=min(depthRange) and f.profile_range<max(depthRange)))
        maxlag = np.ceil(maxOffset/dr)
        [AC.RANGEIND(ii,:),AC.AMPCOR(ii,:),~,AC.LAGS(ii,:)] = fmcw_xcorr(f.specCor,g.specCor,fi,maxlag);
        AC.RANGE(ii,:) = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,AC.RANGEIND(ii,:));
        [~,mci] = max(AC.AMPCOR(ii,:));
        AC.range(ii) = AC.RANGE(ii,mci); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
        AC.lags(ii) = AC.LAGS(ii,mci);
        AC.ampCor(ii) = AC.AMPCOR(ii,mci);
        AC.dh(ii) = dr*AC.lags(ii); % Range offset (m) between segments

        % Quality checks on best correlation
        % Check whether correlation is limited by maxlag
        if mci == 1 || mci == size(AC.LAGS,2) 
            AC.ampCor(ii) = 0;
        end
        % Check prominence of peak (how much better than the next match)
        [cpk,~] = findpeaks(AC.AMPCOR(ii,:),'sortstr','descend','npeaks',2);
        if isempty(cpk) % no peaks!
            AC.ampCorProm(ii) = 0;
        elseif numel(cpk)==1
            AC.ampCorProm(ii) = 1; % this is the only maximum
        else
            AC.ampCorProm(ii) = cpk(1) - cpk(2); % Absolute prominence
        
    
    AC.isGood = AC.ampCor>=cfg.minAmpCor & AC.ampCorProm>=cfg.minAmpCorProm;

    % Now fit a polnomial through the lags
    AC.P = polyfit(AC.range(AC.isGood),AC.lags(AC.isGood),cfg.polyOrder);

def xcorr(f,g,fi,maxlag,p,fe=None,ge=None):
    if len(maxlag)==1:
        lags = np.arange(-abs(maxlag),abs(maxlag)+1,1)
    elif len(maxlag)==2:
        lags = np.arange(maxlag[0],maxlag[1]+1,1)
    
    n = max(abs(lags))
    fi = fi+n;
    zp = zeros(n,1);
    
    
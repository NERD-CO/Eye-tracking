%% First derivative mask

    % samps = samples.copy(deep=True)
    % for col in mask_cols:
    %     d = samps[col].diff()
    %     m = samps[col].diff().mean()
    %     s = samps[col].diff().std() * threshold
    %     # TODO: check this works properly
    %     samps[col] = samps[col].mask((d < (m - s)) | (d > (m + s)))
    %     samps.loc[samps[col] == 0, col] = np.nan
    % return samps
    % 



    %%% Percent change
    %     for col in pupil_cols:
    %     new_col = col + "_pc"
    %     ranges[new_col] = (ranges[col] / baselines[col] - 1).values * 100
    % return ranges



    %%% butterworth
    %     cutoff_freq : float
    %     Normalised cut-off frequency in hertz. For 4 Hz cut-off, this should
    %     4/(sample_rate/2). The default is .01.
    % 
    %         samps = samples if inplace else samples.copy(deep=True)
    % B, A = signal.butter(filt_order, cutoff_freq, output="BA")
    % samps[fields] = samps[fields].apply(
    %     lambda x: signal.filtfilt(B, A, x), axis=0
    % )
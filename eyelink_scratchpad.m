% Trial 1 Start
% Line edf.FEVENT row 74
% Type = 7 for start
% Type = 8 for end

% Trial 1 / Fix 1
stT = edf3.Samples.time(13);
etT = edf3.Samples.time(102);

aGxE1 = mean(edf3.Samples.gx(13:102,1)); 
aGxE2 = mean(edf3.Samples.gx(13:102,2));

aGyE1 = mean(edf3.Samples.gy(13:102,1)); 
aGyE2 = mean(edf3.Samples.gy(13:102,2));

% time
fs = 1000;
msecs = (length(13:102)/fs)*1000;

% pupil data
aPaE1 = mean(edf3.Samples.pa(13:102,1));
aPaE2 = mean(edf3.Samples.pa(13:102,2));
% angle

% edf3.Events.Sfix.time [every other row will give start and stop times for
% edf3.Events.Efix.time 
% fixations

% edf3.Events.Esacc [saccade ] 

% TTLs and I'm set
edf3.RawEdf.FEVENT.message;

% timestamp and convert to number
% 3 column table
% Message, number, timestamp
[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] =...
           Nlx2MatEV('Events.nev', [1 1 1 1 1], 1, 1, [] );

%%

evst2 = EventStrings(contains(EventStrings,'AcqSystem'))
evst3 = extractBetween(evst2,49,54);
evst4 = hex2dec(evst3) 
tabulate(evst4);


%% 97 598
%% 600 1102

%%

find(contains(EventStrings,'0x0000'))

%%
evst2 = EventStrings(contains(EventStrings,'AcqSystem'))
evst3 = extractBetween(evst2,49,54);
evst4 = hex2dec(evst3) 
%%
startLOCs = find(evst4 == 55)


%%

evst2 = EventStrings(contains(EventStrings,'AcqSystem'));
evst3 = extractBetween(evst2,49,54);
evst4 = hex2dec(evst3) 
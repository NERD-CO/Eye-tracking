[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] =...
           Nlx2MatEV('Events.nev', [1 1 1 1 1], 1, 1, [] )

%%
evst22 = EventStrings(1:1011)
evst2 = evst22(contains(evst22,'AcqSystem'))
evst3 = extractBetween(evst2,49,54);
evst4 = hex2dec(evst3) 
tabulate(evst4);

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
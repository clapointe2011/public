bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));

bool monitorConservation(pimple.dict().lookupOrDefault<Switch>("monitorConservation", true));

scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));

scalar startRunIndex(pimple.dict().lookupOrDefault<scalar>("startRunIndex",0));

bool protectOutlet(pimple.dict().lookupOrDefault<Switch>("protectOutlet",false));

bool protectInlet(pimple.dict().lookupOrDefault<Switch>("protectInlet",false));

bool useMULES(pimple.dict().lookupOrDefault<Switch>("useMULES",false));

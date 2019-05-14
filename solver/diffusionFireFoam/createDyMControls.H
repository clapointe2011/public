bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));

scalar nCorrPhi(pimple.dict().lookupOrDefault<scalar>("nCorrPhi",1));

scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));

scalar startRunIndex(pimple.dict().lookupOrDefault<scalar>("startRunIndex",0));

bool protectOutlet(pimple.dict().lookupOrDefault<Switch>("protectOutlet",false));

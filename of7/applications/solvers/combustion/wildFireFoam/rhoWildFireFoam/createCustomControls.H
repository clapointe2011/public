bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));

scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));

scalar startRunIndex(pimple.dict().lookupOrDefault<scalar>("startRunIndex",0));

bool protectMapped(pimple.dict().lookupOrDefault<Switch>("protectMapped",true));

bool addOnlyNearMapped(pimple.dict().lookupOrDefault<Switch>("addOnlyNearMapped",false));

bool useMULES(pimple.dict().lookupOrDefault<Switch>("useMULES",false));

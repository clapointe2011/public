bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));

scalar nCorrPhi(pimple.dict().lookupOrDefault<scalar>("nCorrPhi",1));

scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));

scalar startRunIndex(pimple.dict().lookupOrDefault<scalar>("startRunIndex",0));

//bool protectOutlet(pimple.dict().lookupOrDefault<Switch>("protectOutlet",false));

bool protectMapped(pimple.dict().lookupOrDefault<Switch>("protectMapped",true));

bool addOnlyNearMapped(pimple.dict().lookupOrDefault<Switch>("addOnlyNearMapped",false));

//bool localBlendedU(pimple.dict().lookupOrDefault<Switch>("localBlendedU",false));



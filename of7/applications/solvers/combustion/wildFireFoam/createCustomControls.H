bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));

scalar nCorrPhi(pimple.dict().lookupOrDefault<scalar>("nCorrPhi",1));

scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));

scalar startRunIndex(pimple.dict().lookupOrDefault<scalar>("startRunIndex",0));

bool protectMapped(pimple.dict().lookupOrDefault<Switch>("protectMapped",true));

bool addOnlyNearMapped(pimple.dict().lookupOrDefault<Switch>("addOnlyNearMapped",false));

bool origHEqn(pimple.dict().lookupOrDefault<Switch>("origHEqn", true));

bool origPIMPLE(pimple.dict().lookupOrDefault<Switch>("origPIMPLE", true));


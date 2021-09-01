word ddtSchemeName(U.mesh().ddtScheme(U.name()));
bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));
scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));
bool protectMapped(pimple.dict().lookupOrDefault<Switch>("protectMapped",true));
bool addOnlyNearMapped(pimple.dict().lookupOrDefault<Switch>("addOnlyNearMapped",false));
bool energyInPiso(pimple.dict().lookupOrDefault<Switch>("energyInPiso",false));
bool turbInPiso(pimple.dict().lookupOrDefault<Switch>("turbInPiso",false));
bool combInPimple(pimple.dict().lookupOrDefault<Switch>("combInPimple",true));
bool radInPimple(pimple.dict().lookupOrDefault<Switch>("radInPimple",true));


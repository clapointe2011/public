bool correctPhi(pimple.dict().lookupOrDefault<Switch>("correctPhi", true));

//bool checkMeshCourantNo(pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false));

//bool moveMeshOuterCorrectors(pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false));

scalar startDyMIndex(pimple.dict().lookupOrDefault<scalar>("startDyMIndex",0));

scalar startRunIndex(pimple.dict().lookupOrDefault<scalar>("startRunIndex",0));

bool protectOutlet(pimple.dict().lookupOrDefault<Switch>("protectOutlet",false));

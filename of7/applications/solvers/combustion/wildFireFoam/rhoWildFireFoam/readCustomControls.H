#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

protectMapped = pimple.dict().lookupOrDefault<Switch>("protectMapped",true);

addOnlyNearMapped = pimple.dict().lookupOrDefault<Switch>("addOnlyNearMapped",false);


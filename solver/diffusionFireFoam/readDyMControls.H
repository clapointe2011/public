#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

bool protectOutlet(pimple.dict().lookupOrDefault<Switch>("protectOutlet",false));


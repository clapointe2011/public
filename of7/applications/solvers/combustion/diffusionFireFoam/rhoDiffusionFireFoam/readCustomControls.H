#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

protectOutlet = pimple.dict().lookupOrDefault<Switch>("protectOutlet",false);

monitorConservation = pimple.dict().lookupOrDefault<Switch>("monitorConservation", true);

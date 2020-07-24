#include "fvCFD.H"

class CostFunction{
public:
    CostFunction(Foam::fvMesh& mesh);
    Foam::scalar eval();
private:
    Foam::fvMesh& mesh;
    Foam::scalar eval_incompPressureDrop(const Foam::wordList& patches);
    Foam::scalar eval_flowOutlet(const Foam::wordList& patches);
    Foam::scalar eval_sonicTotalPressure(const Foam::wordList& patches);
    Foam::scalar eval_rhoCentral(const Foam::wordList& patches);
    Foam::scalar eval_reacSonicTotalPressure(const Foam::wordList& patches);
    Foam::scalar eval_flowSection(const Foam::wordList& patches);
    Foam::scalar eval_sonicFlowPoint(const Foam::wordList& patches);
    Foam::scalar eval_reacFlow(const Foam::wordList& patches);
    Foam::scalar eval_rollerTemp(const Foam::wordList& patches);
    Foam::scalar eval_shockAngle(const Foam::wordList& patches);
    Foam::scalar eval_vorticity(const Foam::wordList& patches);
};

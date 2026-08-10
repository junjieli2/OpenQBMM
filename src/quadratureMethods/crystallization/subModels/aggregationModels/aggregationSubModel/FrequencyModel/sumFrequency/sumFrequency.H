#ifndef sumFrequency_H
#define sumFrequency_H

#include "crystalCollisionFrequency.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalCollisionFrequencies
{

class sumFrequency
:
    public crystalCollisionFrequency
{
    PtrList<dictionary> componentDicts_;
    PtrList<crystalCollisionFrequency> components_;

public:

    TypeName("sum");

    sumFrequency
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~sumFrequency();

    virtual void update();

    virtual scalar beta
    (
        const scalar& d1,
        const scalar& d2,
        const vector& Ur,
        const label celli
    ) const;
};

} // End namespace crystalCollisionFrequencies
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

#endif

// ************************************************************************* //

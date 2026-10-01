// Shear-induced aggregate breakage: Xiao et al., CEJ 538 (2026), Eq. 32-33.
#ifndef shearPowerLawBreakup_H
#define shearPowerLawBreakup_H

#include "breakupKernel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
class shearPowerLawBreakup : public breakupKernel
{
    scalar sizeExponent_;
    scalar shearExponent_;
    dimensionedScalar bStar_;
    const volScalarField& nu_;
    const volScalarField& epsilon_;

public:
    TypeName("shearPowerLawBreakup");
    shearPowerLawBreakup(const dictionary& dict, const fvMesh& mesh);
    virtual scalar Kb
    (
        const scalar& abscissa,
        const label celli,
        const label environment = 0
    ) const;
};
}
}
}
#endif

#include "shearPowerLawBreakup.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
defineTypeNameAndDebug(shearPowerLawBreakup, 0);
addToRunTimeSelectionTable(breakupKernel, shearPowerLawBreakup, dictionary);

shearPowerLawBreakup::shearPowerLawBreakup
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    breakupKernel(dict, mesh),
    sizeExponent_(dict.lookupOrDefault<scalar>("sizeExponent", 2.0)),
    shearExponent_(dict.lookupOrDefault<scalar>("shearExponent", 1.85)),
    bStar_
    (
        "bStar",
        dimensionSet(0, -sizeExponent_, shearExponent_ - 1, 0, 0, 0, 0),
        dict
    ),
    nu_(mesh.lookupObject<volScalarField>("nu")),
    epsilon_(mesh.lookupObject<volScalarField>("epsilon"))
{
    if
    (
        !std::isfinite(sizeExponent_) || sizeExponent_ < 0
     || !std::isfinite(shearExponent_) || shearExponent_ < 0
     || !std::isfinite(bStar_.value()) || bStar_.value() < 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "bStar, sizeExponent and shearExponent must be finite "
            << "and non-negative." << exit(FatalIOError);
    }
    if
    (
        nu_.dimensions() != dimViscosity
     || epsilon_.dimensions() != dimVelocity*dimVelocity/dimTime
    )
    {
        FatalIOErrorInFunction(dict)
            << "nu must have dimensions m2/s and epsilon m2/s3."
            << exit(FatalIOError);
    }
}

scalar shearPowerLawBreakup::Kb
(
    const scalar& length,
    const label celli,
    const label environment
) const
{
    const scalar nu = nu_[celli];
    const scalar eps = epsilon_[celli];
    if
    (
        !std::isfinite(nu) || nu <= 0
     || !std::isfinite(eps) || eps < 0
     || !std::isfinite(length) || length < 0
    )
    {
        FatalErrorInFunction
            << "Invalid breakup input in cell " << celli
            << ": nu=" << nu << ", epsilon=" << eps
            << ", L=" << length << exit(FatalError);
    }
    // Explicit exponent-zero limits also permit constant-rate regression tests.
    const scalar shearFactor = shearExponent_ == 0
        ? 1.0 : Foam::pow(Foam::sqrt(eps/nu), shearExponent_);
    const scalar sizeFactor = sizeExponent_ == 0
        ? 1.0 : Foam::pow(length, sizeExponent_);
    const scalar rate = bStar_.value()*shearFactor*sizeFactor;
    if (!std::isfinite(rate))
    {
        FatalErrorInFunction
            << "Non-finite breakup rate in cell " << celli
            << exit(FatalError);
    }
    return rate;
}
}
}
}

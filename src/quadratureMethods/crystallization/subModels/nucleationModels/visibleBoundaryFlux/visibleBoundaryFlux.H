/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is distributed under the GNU General Public License v3 or later.

Class
    Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux

Description
    Transfers an unresolved log-normal crystal population into the resolved
    (optically visible) population. The prescribed number flux is a gamma
    pulse

        J(t) = A*x^p*exp(-x), x = (t - startTime)/timeScale,

    for x >= 0 and an integer p >= 0. This is a population-transfer source,
    not physical nucleation: its M3 source must not be removed from dissolved
    solute a second time. The model writes an unresolved-M3 reservoir whose
    decrease matches the cumulative resolved M3 added by this source.

SourceFiles
    visibleBoundaryFlux.C

\*---------------------------------------------------------------------------*/

#ifndef visibleBoundaryFlux_H
#define visibleBoundaryFlux_H

#include "nucleationModel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

class visibleBoundaryFlux
:
    public nucleationModel
{
    // Private data

        //- Peak-scale coefficient of the number flux [1/(m3 s)]
        dimensionedScalar amplitude_;

        //- Absolute simulation time at which the flux starts [s]
        scalar startTime_;

        //- Gamma-pulse time scale [s]
        scalar timeScale_;

        //- Non-negative integer gamma-pulse exponent
        label shapeExponent_;

        //- Median entry size [m]
        dimensionedScalar inletMedian_;

        //- Log-space standard deviation of the entry population
        scalar inletSigmaLog_;

        //- Initial unresolved third moment (dimensionless)
        scalar initialUnresolvedM3_;

        //- Instantaneous unresolved-to-visible number flux
        mutable volScalarField fluxRate_;

        //- Remaining unresolved third moment
        mutable volScalarField unresolvedM3_;

        //- Cumulative M3 transferred into the visible population
        mutable volScalarField cumulativeTransferredM3_;

        //- Last global time index used to update the discrete reservoir
        mutable label lastTimeIndex_;

        //- Discretely accumulated transferred M3
        mutable scalar cumulativeTransferredM3Value_;


    // Private member functions

        scalar factorial(const label n) const;

        scalar inletRawMoment(const label order) const;

        scalar numberFlux() const;

        scalar cumulativeTransferredNumber() const;

        void updateDiagnostics(const label celli) const;


public:

    TypeName("visibleBoundaryFlux");

    visibleBoundaryFlux
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~visibleBoundaryFlux();

    //- This transfers an existing crystal population and therefore must not
    //  contribute to the dissolved-solute source.
    virtual bool contributesToSpeciesTransfer() const
    {
        return false;
    }

    virtual scalar nucleationSource
    (
        const label& momentOrder,
        const label celli,
        const label environment = 0
    ) const;
};

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

#endif

// ************************************************************************* //

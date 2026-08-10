/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenQBMM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenQBMM and is distributed under the GNU General
    Public License, version 3 or later.
\*---------------------------------------------------------------------------*/

#include "borsosPrimaryNucleation.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(borsosPrimaryNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        borsosPrimaryNucleation,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::nucleationModels::
borsosPrimaryNucleation::borsosPrimaryNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    kp0_("kp0", inv(dimVolume*dimTime), dict),
    Ep_("Ep", dimEnergy/dimMoles, dict),
    R_
    (
        dict.lookupOrDefault
        (
            "R",
            dimensionedScalar
            (
                "R",
                dimEnergy/dimMoles/dimTemperature,
                8.31446261815324
            )
        )
    ),
    ke_(readScalar(dict.lookup("ke"))),
    nucleationScale_(dict.lookupOrDefault<scalar>("nucleationScale", 1.0)),
    T_
    (
        mesh.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("TField", "T")
        )
    ),
    C_
    (
        mesh.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("soluteField", "C")
        )
    ),
    Csat_
    (
        mesh.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("saturationField", "Csat")
        )
    ),
    rate_
    (
        IOobject
        (
            dict.lookupOrDefault<word>
            (
                "rateField",
                "borsosNucleationRate"
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", inv(dimVolume*dimTime), 0.0)
    ),
    Lnucleation_
    (
        "Lnucleation",
        dimLength,
        dict.lookupOrDefault<scalar>
        (
            "Lnucleation",
            dict.lookupOrDefault<scalar>("d_nucleation", 1.0e-6)
        )
    )
{
    if
    (
        kp0_.value() < 0
     || Ep_.value() < 0
     || R_.value() <= 0
     || ke_ < 0
     || nucleationScale_ < 0
     || Lnucleation_.value() <= 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "kp0, Ep, ke and nucleationScale must be non-negative; "
            << "R and Lnucleation must be positive." << exit(FatalIOError);
    }
}


Foam::populationBalanceSubModels::nucleationModels::
borsosPrimaryNucleation::~borsosPrimaryNucleation()
{}


Foam::scalar Foam::populationBalanceSubModels::nucleationModels::
borsosPrimaryNucleation::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    (void)environment;

    scalar Bp = 0.0;
    const scalar temperature = T_[celli];
    const scalar csat = Csat_[celli];
    const scalar concentration = C_[celli];

    if (temperature > 0 && csat > SMALL && concentration > csat)
    {
        const scalar logS = std::log(concentration/csat);
        if (logS > SMALL)
        {
            Bp =
                nucleationScale_*kp0_.value()
               *std::exp(-Ep_.value()/(R_.value()*temperature))
               *std::exp(-ke_/sqr(logS));
        }
    }

    if (!std::isfinite(Bp) || Bp < 0)
    {
        FatalErrorInFunction
            << "Borsos primary nucleation produced invalid Bp=" << Bp
            << " for T=" << temperature << ", C=" << concentration
            << " and Csat=" << csat << exit(FatalError);
    }

    const_cast<volScalarField&>(rate_)[celli] = Bp;

    return Bp*std::pow(Lnucleation_.value(), momentOrder);
}

// ************************************************************************* //

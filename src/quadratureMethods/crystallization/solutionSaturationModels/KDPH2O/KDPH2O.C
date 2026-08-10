/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM and is distributed under the
    GNU General Public License, version 3 or later.

\*---------------------------------------------------------------------------*/

#include "KDPH2O.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

namespace Foam
{
namespace solutionSaturationModels
{
    defineTypeNameAndDebug(KDPH2O, 0);
    addToRunTimeSelectionTable(solutionSaturationModel, KDPH2O, dictionary);
}
}


Foam::solutionSaturationModels::KDPH2O::KDPH2O
(
    const dictionary& dict,
    const objectRegistry& db
)
:
    solutionSaturationModel(db),
    rhoSolvent_
    (
        dimensionedScalar::getOrDefault
        (
            "rhoSolvent",
            dict,
            dimDensity,
            1000.0
        )
    )
{
    if
    (
        !std::isfinite(rhoSolvent_.value())
     || rhoSolvent_.value() <= 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "rhoSolvent must be finite and positive, but is "
            << rhoSolvent_ << exit(FatalIOError);
    }
}


Foam::scalar
Foam::solutionSaturationModels::KDPH2O::saturationConcentration
(
    const scalar temperatureKelvin
) const
{
    if (!std::isfinite(temperatureKelvin) || temperatureKelvin <= 0)
    {
        FatalErrorInFunction
            << "KDPH2O requires a finite positive temperature in K, but "
            << "received " << temperatureKelvin << exit(FatalError);
    }

    const scalar temperatureCelsius = temperatureKelvin - 273.15;
    const scalar massRatio =
        9.3027e-5*sqr(temperatureCelsius)
      - 9.7629e-5*temperatureCelsius
      + 0.2087;
    const scalar concentration = rhoSolvent_.value()*massRatio;

    if (!std::isfinite(concentration) || concentration <= 0)
    {
        FatalErrorInFunction
            << "KDPH2O produced an invalid saturation concentration "
            << concentration << " kg/m3 at T = " << temperatureKelvin
            << " K. The published correlation is intended for aqueous KDP "
            << "over its experimental temperature range."
            << exit(FatalError);
    }

    return concentration;
}


Foam::tmp<Foam::volScalarField>
Foam::solutionSaturationModels::KDPH2O::Csat
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tCsat
    (
        volScalarField::New
        (
            "Csat",
            T.mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );

    volScalarField& Csat = tCsat.ref();

    forAll(Csat, celli)
    {
        Csat[celli] = saturationConcentration(T[celli]);
    }

    volScalarField::Boundary& CsatBf = Csat.boundaryFieldRef();

    forAll(CsatBf, patchi)
    {
        scalarField& CsatPatch = CsatBf[patchi];
        const scalarField& temperaturePatch = T.boundaryField()[patchi];

        forAll(CsatPatch, facei)
        {
            CsatPatch[facei] =
                saturationConcentration(temperaturePatch[facei]);
        }
    }

    return tCsat;
}


// ************************************************************************* //

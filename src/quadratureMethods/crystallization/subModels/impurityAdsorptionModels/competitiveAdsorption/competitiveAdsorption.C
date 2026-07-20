/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "competitiveAdsorption.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace impurityAdsorptionModels
{
    defineTypeNameAndDebug(competitiveAdsorption, 0);

    addToRunTimeSelectionTable
    (
        impurityAdsorptionModel,
        competitiveAdsorption,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::impurityAdsorptionModels
::competitiveAdsorption::competitiveAdsorption
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    impurityAdsorptionModel(dict, mesh),
    Ki_("Ki", inv(dimDensity), dict),
    Ks_("Ks", inv(dimDensity), dict)
{
    const dimensionSet inverseDensity(inv(dimDensity));

    if (Ki_.dimensions() != inverseDensity)
    {
        FatalErrorInFunction
            << "Ki has dimensions " << Ki_.dimensions()
            << ", expected " << inverseDensity
            << exit(FatalError);
    }

    if (Ks_.dimensions() != inverseDensity)
    {
        FatalErrorInFunction
            << "Ks has dimensions " << Ks_.dimensions()
            << ", expected " << inverseDensity
            << exit(FatalError);
    }

    if (!std::isfinite(Ki_.value()) || Ki_.value() < 0)
    {
        FatalErrorInFunction
            << "Ki must be finite and non-negative, but is " << Ki_.value()
            << exit(FatalError);
    }

    if (!std::isfinite(Ks_.value()) || Ks_.value() < 0)
    {
        FatalErrorInFunction
            << "Ks must be finite and non-negative, but is " << Ks_.value()
            << exit(FatalError);
    }

    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::impurityAdsorptionModels
::competitiveAdsorption::adsorbedFraction
(
    const scalar soluteConcentration,
    const scalar impurityConcentration
) const
{
    if
    (
        !std::isfinite(soluteConcentration)
     || !std::isfinite(impurityConcentration)
    )
    {
        FatalErrorInFunction
            << "Competitive adsorption requires finite concentrations. "
            << "Received C = " << soluteConcentration
            << " and Ci = " << impurityConcentration
            << exit(FatalError);
    }

    // Accommodate small transport-solver undershoots without hiding a
    // physically invalid concentration field.
    const scalar soluteTolerance =
        ROOTSMALL*max(scalar(1), mag(soluteConcentration));
    const scalar impurityTolerance =
        ROOTSMALL*max(scalar(1), mag(impurityConcentration));

    if (soluteConcentration < -soluteTolerance)
    {
        FatalErrorInFunction
            << "Solute concentration C must be non-negative. Received "
            << soluteConcentration << ", below the round-off tolerance "
            << soluteTolerance
            << exit(FatalError);
    }

    if (impurityConcentration < -impurityTolerance)
    {
        FatalErrorInFunction
            << "Impurity concentration Ci must be non-negative. Received "
            << impurityConcentration << ", below the round-off tolerance "
            << impurityTolerance
            << exit(FatalError);
    }

    const scalar boundedSolute = max(soluteConcentration, scalar(0));
    const scalar boundedImpurity = max(impurityConcentration, scalar(0));
    const scalar numerator = Ki_.value()*boundedImpurity;
    const scalar denominator =
        1.0 + numerator + Ks_.value()*boundedSolute;

    if (!std::isfinite(numerator) || !std::isfinite(denominator))
    {
        FatalErrorInFunction
            << "Competitive adsorption overflowed for C = "
            << boundedSolute << " and Ci = " << boundedImpurity
            << ": numerator = " << numerator
            << ", denominator = " << denominator
            << exit(FatalError);
    }

    if (!(denominator > VSMALL))
    {
        FatalErrorInFunction
            << "Competitive adsorption denominator must be positive. "
            << "Calculated " << denominator << " from C = "
            << soluteConcentration << " and Ci = "
            << impurityConcentration
            << exit(FatalError);
    }

    const scalar theta = numerator/denominator;

    if (!std::isfinite(theta))
    {
        FatalErrorInFunction
            << "Competitive adsorption produced non-finite theta = "
            << theta << " for numerator = " << numerator
            << " and denominator = " << denominator
            << exit(FatalError);
    }

    if (theta < -ROOTSMALL || theta > 1.0 + ROOTSMALL)
    {
        FatalErrorInFunction
            << "Competitive adsorption produced theta = " << theta
            << ", outside [0, 1]. Check that C and Ci are non-negative."
            << exit(FatalError);
    }

    // Remove only round-off excursions at the physical bounds.
    return min(max(theta, scalar(0)), scalar(1));
}


void Foam::populationBalanceSubModels::impurityAdsorptionModels
::competitiveAdsorption::correct()
{
    forAll(thetaImpurity_, celli)
    {
        thetaImpurity_[celli] =
            adsorbedFraction(solute_[celli], impurity_[celli]);
    }

    volScalarField::Boundary& thetaBf =
        thetaImpurity_.boundaryFieldRef();

    forAll(thetaBf, patchi)
    {
        scalarField& thetaPatch = thetaBf[patchi];
        const scalarField& solutePatch = solute_.boundaryField()[patchi];
        const scalarField& impurityPatch = impurity_.boundaryField()[patchi];

        forAll(thetaPatch, facei)
        {
            thetaPatch[facei] =
                adsorbedFraction
                (
                    solutePatch[facei],
                    impurityPatch[facei]
                );
        }
    }
}


// ************************************************************************* //

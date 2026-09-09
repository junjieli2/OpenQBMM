/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Junjie Li
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HounslowEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalAggregationEfficiencies
{
    defineTypeNameAndDebug(Hounslow, 0);

    addToRunTimeSelectionTable
    (
        crystalAggregationEfficiency,
        Hounslow,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::Hounslow
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crystalAggregationEfficiency(dict, mesh),
    Mg_(dict.lookupOrDefault<scalar>("Mg", 1.16)),
    sigmag_(dict.lookupOrDefault<scalar>("sigmag", 3.08)),
    Lstar_sigmaY_(dict.lookupOrDefault<scalar>("Lstar_sigmaY", 1.35)),
    rhoFluid_(dict.lookupOrDefault<scalar>("rhoFluid", 1000.0)),
    qBlend_(max(dict.lookupOrDefault<scalar>("qBlend", 0.05), SMALL)),
    qMin_(max(dict.lookupOrDefault<scalar>("qMin", 0.05), SMALL)),
    useSizeRatioQ_(dict.lookupOrDefault<bool>("useSizeRatioQ", true)),
    usePairGrowth_(dict.lookupOrDefault<bool>("usePairGrowth", false)),
    pairCgGrowth_
    (
        dict.lookupOrDefault
        (
            "pairCgGrowth",
            dimensionedScalar("pairCgGrowth", dimLength/dimTime, 1e-6)
        )
    ),
    pairGrowthPow_(dict.lookupOrDefault<scalar>("pairGrowthPow", 1.0)),
    pairSigmaGrowthCrit_
    (
        dict.lookupOrDefault<scalar>("pairSigmaGrowthCrit", 0.0)
    ),
    pairSigmaMax_(dict.lookupOrDefault<scalar>("pairSigmaMax", 10.0)),
    pairMaxGrowth_
    (
        dict.lookupOrDefault
        (
            "pairMaxGrowth",
            dimensionedScalar("pairMaxGrowth", dimLength/dimTime, 1e-6)
        )
    ),
    pairLc_
    (
        dict.lookupOrDefault
        (
            "pairLc",
            dimensionedScalar("pairLc", dimLength, 200e-6)
        )
    ),
    pairQSize_(dict.lookupOrDefault<scalar>("pairQSize", 2.0)),
    useDimensionlessStrength_
    (
        dict.lookupOrDefault<bool>("useDimensionlessStrength", false)
    ),
    PiC_(dict.lookupOrDefault<scalar>("PiC", 1.0)),
    referenceLength_
    (
        max(dict.lookupOrDefault<scalar>("referenceLength", 1e-4), SMALL)
    ),
    referenceGrowthRate_
    (
        max(dict.lookupOrDefault<scalar>("referenceGrowthRate", 1e-6), SMALL)
    ),
    referenceShearRate_
    (
        max(dict.lookupOrDefault<scalar>("referenceShearRate", 290.4737509655563), SMALL)
    ),
    referenceViscosity_
    (
        max(dict.lookupOrDefault<scalar>("referenceViscosity", 1e-3), SMALL)
    ),
    GField_
    (
        mesh.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("growthRateField", "growthRate")
        )
    ),
    sigma_
    (
        mesh.lookupObject<volScalarField>("sigma")
    ),
    nu_
    (
        mesh.lookupObject<volScalarField>("nu")
    ),
    epsilon_
    (
        mesh.lookupObject<volScalarField>("epsilon")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::~Hounslow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::pairGrowthRate
(
    const scalar contactLength,
    const label celli
) const
{
    const scalar localSigma = sigma_[celli];
    if (!std::isfinite(localSigma))
    {
        return 0.0;
    }

    const scalar sigmaEff =
        max(min(localSigma, pairSigmaMax_) - pairSigmaGrowthCrit_, 0.0);
    if (sigmaEff <= SMALL)
    {
        return 0.0;
    }

    const scalar L = max(contactLength, SMALL);
    const scalar Gbase =
        pairCgGrowth_.value()*Foam::pow(sigmaEff, pairGrowthPow_);
    scalar sizeFactor =
        1.0/(1.0 + Foam::pow(L/max(pairLc_.value(), SMALL), pairQSize_));
    sizeFactor = max(min(sizeFactor, 1.0), 0.0);

    const scalar G = min(Gbase*sizeFactor, pairMaxGrowth_.value());
    return std::isfinite(G) && G > SMALL ? G : scalar(0.0);
}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
crystalAggregationEfficiencies::Hounslow::Pc
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    (void)Ur;

    if (!std::isfinite(d1) || !std::isfinite(d2) || d1 <= SMALL || d2 <= SMALL)
    {
        return 0.0;
    }

    const scalar nu = max(nu_[celli], SMALL);
    const scalar eps = max(epsilon_[celli], SMALL);
    if (!std::isfinite(nu) || !std::isfinite(eps))
    {
        return 0.0;
    }

    const scalar gammaDot = Foam::sqrt(eps/nu);
    if (!std::isfinite(gammaDot) || gammaDot <= SMALL)
    {
        return 0.0;
    }

    const scalar lambda = max(d1, d2)/max(min(d1, d2), SMALL);
    if (!std::isfinite(lambda) || lambda <= SMALL)
    {
        return 0.0;
    }

    const scalar q = useSizeRatioQ_
      ? 1.0 - 0.328*mag(Foam::log(lambda))
      : 1.0;
    if (!std::isfinite(q) || q <= 0.0)
    {
        return 0.0;
    }

    const scalar blendFraction = min(1.0, max(0.0, q/qBlend_));
    const scalar qTaper = sqr(blendFraction)*(3.0 - 2.0*blendFraction);
    const scalar L = max(q, qMin_)*Foam::sqrt(d1*d2);
    if (!std::isfinite(L) || L <= SMALL || !std::isfinite(qTaper))
    {
        return 0.0;
    }

    const scalar G = usePairGrowth_
      ? pairGrowthRate(L, celli)
      : max(GField_[celli], scalar(0.0));
    if (!std::isfinite(G) || G <= SMALL)
    {
        return 0.0;
    }

    const scalar mu = rhoFluid_*nu;
    if (!std::isfinite(mu) || mu <= SMALL || sigmag_ <= 1.0)
    {
        return 0.0;
    }

    scalar strength = 0.0;
    if (useDimensionlessStrength_)
    {
        const scalar shearRatio = gammaDot/referenceShearRate_;
        const scalar lengthRatio = L/referenceLength_;
        const scalar denominator =
            sqr(shearRatio)*(mu/referenceViscosity_)*sqr(lengthRatio);
        strength = PiC_*(G/referenceGrowthRate_)/(denominator + SMALL);
    }
    else
    {
        const scalar M =
            (Lstar_sigmaY_*G)/(sqr(gammaDot)*mu*sqr(L) + SMALL);
        strength = Mg_ > SMALL ? M/Mg_ : scalar(0.0);
    }

    if (!std::isfinite(strength) || strength <= SMALL)
    {
        return 0.0;
    }

    const scalar arg =
        Foam::log(strength)/(Foam::sqrt(2.0)*Foam::log(sigmag_));
    if (!std::isfinite(arg))
    {
        return 0.0;
    }

    const scalar psi = qTaper*0.5*(1.0 + Foam::erf(arg));
    if (!std::isfinite(psi))
    {
        return 0.0;
    }

    return max(0.0, min(1.0, psi));
}


// ************************************************************************* //

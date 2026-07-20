/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
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

#include "areaNucleation.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace
{

Foam::word resolveAreaMomentName
(
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh,
    const Foam::word& momentGroup
)
{
    if (dict.found("m2Name"))
    {
        const Foam::word configuredName(dict.lookup("m2Name"));

        if (!mesh.foundObject<Foam::volScalarField>(configuredName))
        {
            FatalIOErrorInFunction(dict)
                << "Configured second-moment field '" << configuredName
                << "' was not found in the mesh object registry."
                << exit(Foam::FatalIOError);
        }

        return configuredName;
    }

    const Foam::word bivariateName
    (
        Foam::IOobject::groupName
        (
            Foam::IOobject::groupName("moment", "20"),
            momentGroup
        )
    );

    if (mesh.foundObject<Foam::volScalarField>(bivariateName))
    {
        return bivariateName;
    }

    const Foam::word legacyName
    (
        Foam::IOobject::groupName
        (
            Foam::IOobject::groupName("moment", "2"),
            momentGroup
        )
    );

    if (mesh.foundObject<Foam::volScalarField>(legacyName))
    {
        return legacyName;
    }

    FatalIOErrorInFunction(dict)
        << "No second-moment field was found. Tried bivariate field '"
        << bivariateName << "' and legacy field '" << legacyName << "'. "
        << "Set m2Name explicitly to use another field."
        << exit(Foam::FatalIOError);

    return bivariateName;
}

}

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(areaNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        areaNucleation,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::nucleationModels::areaNucleation::areaNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    kb_("kb", inv(dimVolume*dimTime), dict),
    b_(readScalar(dict.lookup("b"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    momentGroup_(dict.lookupOrDefault<word>("momentGroup", "populationBalance")),
    m2Name_(resolveAreaMomentName(dict, mesh, momentGroup_)),
    sigma_
    (
        mesh.lookupObject<volScalarField>("sigma")
    ),
    m2_
    (
        mesh.lookupObject<volScalarField>(m2Name_)
    ),
    JField_
    (
        IOobject
        (
            "nucleationRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", inv(dimVolume*dimTime), 0.0)
    ),
    d_nucleation_
    (
        "d_nucleation",
        dimLength,
        dict.lookupOrDefault<scalar>("d_nucleation", dict.lookupOrDefault<scalar>("deltaWidth", 1.0e-6))
    )

{}

Foam::populationBalanceSubModels::nucleationModels::areaNucleation::~areaNucleation()
{}

Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::areaNucleation::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    const scalar sigmaVal = sigma_[celli];
    const scalar d_nucleation = d_nucleation_.value();
    const scalar areaInhibition =
        std::exp(-alpha_ * max(m2_[celli], scalar(0)));

    const scalar J =
        kb_.value()
      * std::pow(max(sigmaVal, scalar(0)), b_)
      * areaInhibition;

    const scalar nucleationSource = J * std::pow(d_nucleation, momentOrder);

    const_cast<volScalarField&>(JField_)[celli] = J;

    return nucleationSource;
}

// ************************************************************************* //

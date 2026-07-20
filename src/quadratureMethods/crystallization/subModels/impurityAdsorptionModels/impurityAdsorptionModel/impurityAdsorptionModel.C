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

#include "impurityAdsorptionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(impurityAdsorptionModel, 0);
    defineRunTimeSelectionTable(impurityAdsorptionModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::impurityAdsorptionModel
::impurityAdsorptionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    soluteFieldName_(dict.lookupOrDefault<word>("soluteField", "C")),
    impurityFieldName_(dict.lookupOrDefault<word>("impurityField", "Ci")),
    solute_(mesh.lookupObject<volScalarField>(soluteFieldName_)),
    impurity_(mesh.lookupObject<volScalarField>(impurityFieldName_)),
    thetaImpurity_
    (
        IOobject
        (
            "thetaImpurity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::calculatedType()
    )
{
    if (solute_.dimensions() != dimDensity)
    {
        FatalErrorInFunction
            << "Solute field '" << soluteFieldName_ << "' has dimensions "
            << solute_.dimensions() << ", expected " << dimDensity
            << exit(FatalError);
    }

    if (impurity_.dimensions() != dimDensity)
    {
        FatalErrorInFunction
            << "Impurity field '" << impurityFieldName_ << "' has dimensions "
            << impurity_.dimensions() << ", expected " << dimDensity
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr
<Foam::populationBalanceSubModels::impurityAdsorptionModel>
Foam::populationBalanceSubModels::impurityAdsorptionModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word modelType
    (
        dict.lookupOrDefault<word>
        (
            "impurityAdsorptionModel",
            "competitiveAdsorption"
        )
    );

    Info<< "Selecting impurityAdsorptionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInLookup
        (
            dict,
            "impurityAdsorptionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<impurityAdsorptionModel>(cstrIter()(dict, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::impurityAdsorptionModel::preUpdate()
{
    correct();
}


// ************************************************************************* //

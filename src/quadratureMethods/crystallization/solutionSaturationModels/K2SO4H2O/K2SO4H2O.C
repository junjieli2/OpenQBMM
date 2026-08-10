/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "K2SO4H2O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{
    defineTypeNameAndDebug(K2SO4H2O, 0);
    addToRunTimeSelectionTable(solutionSaturationModel, K2SO4H2O, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::K2SO4H2O::K2SO4H2O
(
    const dictionary& dict,
    const objectRegistry& db
)
:
    solutionSaturationModel(db),
    A_(dict.lookupOrDefault<scalar>("A", -249.369)),
    B_(dict.lookupOrDefault<scalar>("B", 0.3761)),
    C_(dict.lookupOrDefault<scalar>("C", 0.0029415))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionSaturationModels::K2SO4H2O::~K2SO4H2O()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::solutionSaturationModels::K2SO4H2O::Csat
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

    // 1. Calculate Internal Field
    forAll(Csat, celli)
    {
        scalar TK = T[celli]; 
        
        Csat[celli] = A_ + B_*TK + C_*pow(TK, 2);
    }

    // 2. Calculate Boundary Field
    volScalarField::Boundary& CsatBf = Csat.boundaryFieldRef();

    forAll(Csat.boundaryField(), patchi)
    {
        scalarField& Csatp = CsatBf[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(Csatp, facei)
        {
            scalar TK = Tp[facei];
            Csatp[facei] = A_ + B_*TK + C_*pow(TK, 2);
        }
    }

    return tCsat;
}


// ************************************************************************* //

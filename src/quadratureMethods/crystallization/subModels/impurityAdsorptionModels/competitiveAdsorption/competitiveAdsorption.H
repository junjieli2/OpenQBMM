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

Class
    Foam::populationBalanceSubModels::impurityAdsorptionModels::competitiveAdsorption

Description
    Competitive Langmuir adsorption model

        theta = Ki*Ci/(1 + Ki*Ci + Ks*C)

\*---------------------------------------------------------------------------*/

#ifndef competitiveAdsorption_H
#define competitiveAdsorption_H

#include "impurityAdsorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace impurityAdsorptionModels
{

/*---------------------------------------------------------------------------*\
                Class competitiveAdsorption Declaration
\*---------------------------------------------------------------------------*/

class competitiveAdsorption
:
    public impurityAdsorptionModel
{
    // Private data

        //- Impurity adsorption equilibrium coefficient
        dimensionedScalar Ki_;

        //- Solute adsorption equilibrium coefficient
        dimensionedScalar Ks_;


    // Private Member Functions

        //- Evaluate and validate the local algebraic adsorption fraction
        scalar adsorbedFraction
        (
            const scalar soluteConcentration,
            const scalar impurityConcentration
        ) const;


public:

    //- Runtime type information
    TypeName("competitiveAdsorption");


    // Constructors

        //- Construct from a coefficient dictionary and mesh
        competitiveAdsorption
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~competitiveAdsorption() = default;


    // Member Functions

        //- Recalculate competitive adsorption equilibrium
        virtual void correct() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace impurityAdsorptionModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

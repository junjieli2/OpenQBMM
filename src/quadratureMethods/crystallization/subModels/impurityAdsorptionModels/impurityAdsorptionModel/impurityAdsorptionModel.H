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
    Foam::populationBalanceSubModels::impurityAdsorptionModel

Description
    Runtime-selectable algebraic impurity adsorption model.

SourceFiles
    impurityAdsorptionModel.C

\*---------------------------------------------------------------------------*/

#ifndef impurityAdsorptionModel_H
#define impurityAdsorptionModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                Class impurityAdsorptionModel Declaration
\*---------------------------------------------------------------------------*/

class impurityAdsorptionModel
{
protected:

    // Protected data

        //- Model coefficients
        const dictionary& dict_;

        //- Computational mesh
        const fvMesh& mesh_;

        //- Dissolved solute field name
        const word soluteFieldName_;

        //- Dissolved impurity field name
        const word impurityFieldName_;

        //- Dissolved solute concentration
        const volScalarField& solute_;

        //- Dissolved impurity concentration
        const volScalarField& impurity_;

        //- Adsorbed surface fraction
        volScalarField thetaImpurity_;


public:

    //- Runtime type information
    TypeName("impurityAdsorptionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        impurityAdsorptionModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from a coefficient dictionary and mesh
        impurityAdsorptionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow copy construction
        impurityAdsorptionModel(const impurityAdsorptionModel&) = delete;


    // Selectors

        //- Select the adsorption model (competitiveAdsorption by default)
        static autoPtr<impurityAdsorptionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~impurityAdsorptionModel() = default;


    // Member Functions

        //- Update model state before population-balance source evaluation
        virtual void preUpdate();

        //- Recalculate the algebraic adsorption state
        virtual void correct() = 0;

        //- Adsorbed fraction in a cell
        scalar theta(const label celli) const
        {
            return thetaImpurity_[celli];
        }

        //- Adsorbed fraction field
        const volScalarField& thetaField() const
        {
            return thetaImpurity_;
        }

        //- Dissolved solute concentration field
        const volScalarField& solute() const
        {
            return solute_;
        }

        //- Dissolved impurity concentration field
        const volScalarField& impurity() const
        {
            return impurity_;
        }


    // Member Operators

        //- Disallow assignment
        void operator=(const impurityAdsorptionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

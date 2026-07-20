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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::basicScalarFieldMomentInversion

Description
    Invert a multivariate scalar moment field using a run-time selectable
    multivariate moment inversion.

SourceFiles
    basicScalarFieldMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef basicScalarFieldMomentInversion_H
#define basicScalarFieldMomentInversion_H

#include "fieldMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multivariateMomentInversion;

/*---------------------------------------------------------------------------*\
             Class basicScalarFieldMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class basicScalarFieldMomentInversion
:
    public fieldMomentInversion
{
    // Private data

        //- Multivariate moment inversion method
        autoPtr<multivariateMomentInversion> momentInverter_;


public:

    //- Runtime type information
    TypeName("basicScalarFieldMomentInversion");


    // Constructors

        //- Construct from dictionary
        basicScalarFieldMomentInversion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );

        //- Disallow default bitwise copy construct
        basicScalarFieldMomentInversion
        (
            const basicScalarFieldMomentInversion&
        ) = delete;


    //- Destructor
    virtual ~basicScalarFieldMomentInversion();


    // Member Functions

        // Scalar moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            );

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            );

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            );


        // Velocity moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes
            );

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes
            );

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            );

        //- Return the smallest accepted value for m0
        virtual scalar smallM0() const;

        //- Return the smallest accepted value for zeta_k
        virtual scalar smallZeta() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicScalarFieldMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::univariateAdvection::multivariateFirstOrderKineticScalar

Description
    First-order kinetic advection for a scalar quadrature with an arbitrary
    number of internal-coordinate dimensions.

SourceFiles
    multivariateFirstOrderKineticScalar.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateFirstOrderKineticScalar_H
#define multivariateFirstOrderKineticScalar_H

#include "univariateMomentAdvection.H"
#include "basicScalarFieldMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace univariateAdvection
{

/*---------------------------------------------------------------------------*\
           Class multivariateFirstOrderKineticScalar Declaration
\*---------------------------------------------------------------------------*/

class multivariateFirstOrderKineticScalar
:
    public univariateMomentAdvection
{
    // Private data

        //- Cell-centred quadrature used for advection
        autoPtr<mappedPtrList<volScalarNode>> nodes_;

        //- Neighbour-side face quadrature
        autoPtr<mappedPtrList<surfaceScalarNode>> nodesNei_;

        //- Owner-side face quadrature
        autoPtr<mappedPtrList<surfaceScalarNode>> nodesOwn_;

        //- Neighbour-side face moments
        surfaceScalarMomentFieldSet momentsNei_;

        //- Owner-side face moments
        surfaceScalarMomentFieldSet momentsOwn_;

        //- Cell moment inversion method
        autoPtr<basicScalarFieldMomentInversion> momentFieldInverter_;

        //- Number of configured quadrature nodes
        label nNodes_;


    // Private Member Functions

        //- Interpolate all weights and scalar abscissae to faces
        void interpolateNodes();


public:

    //- Runtime type information
    TypeName("multivariateFirstOrderKineticScalar");


    // Constructors

        //- Construct from dictionary and scalar quadrature
        multivariateFirstOrderKineticScalar
        (
            const dictionary& dict,
            const scalarQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const supportType& support
        );


    //- Destructor
    virtual ~multivariateFirstOrderKineticScalar();


    // Member Functions

        //- Return the maximum Courant number ensuring realizability
        virtual scalar realizableCo() const;

        //- Update moment advection
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace univariateAdvection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

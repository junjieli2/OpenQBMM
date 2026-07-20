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

\*---------------------------------------------------------------------------*/

#include "basicScalarFieldMomentInversion.H"
#include "addToRunTimeSelectionTable.H"
#include "multivariateMomentInversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicScalarFieldMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        fieldMomentInversion,
        basicScalarFieldMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::basicScalarFieldMomentInversion::basicScalarFieldMomentInversion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
:
    fieldMomentInversion
    (
        dict,
        mesh,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    ),
    momentInverter_
    (
        multivariateMomentInversion::New
        (
            dict.subDict("basicScalarMomentInversion"),
            momentOrders,
            nodeIndexes,
            labelList(1, -1)
        )
    )
{}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::basicScalarFieldMomentInversion::~basicScalarFieldMomentInversion()
{}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::basicScalarFieldMomentInversion::invert
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    const volScalarField& m0(moments(0));

    forAll(m0, celli)
    {
        invertLocalMoments(moments, nodes, celli);
    }

    invertBoundaryMoments(moments, nodes);
}


void Foam::basicScalarFieldMomentInversion::invertBoundaryMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    const volScalarField::Boundary& bf = moments[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            multivariateMomentSet momentsToInvert
            (
                moments.size(),
                momentOrders_,
                moments.supports(),
                smallM0(),
                smallZeta()
            );

            forAll(momentOrders_, momenti)
            {
                const labelList& momentOrder = momentOrders_[momenti];

                momentsToInvert(momentOrder) =
                    moments(momentOrder).boundaryField()[patchi][facei];
            }

            if (!momentInverter_().invert(momentsToInvert))
            {
                FatalErrorInFunction
                    << "Moment set on patch " << bf[patchi].patch().name()
                    << ", face " << facei << " is not realizable"
                    << abort(FatalError);
            }

            const mappedScalarList& weights(momentInverter_().weights());
            const mappedList<scalarList>& abscissae
            (
                momentInverter_().abscissae()
            );

            forAll(nodes, nodei)
            {
                const labelList& nodeIndex = nodeIndexes_[nodei];
                volScalarNode& node = nodes[nodei];

                node.weight().boundaryFieldRef()[patchi][facei] =
                    weights(nodeIndex);

                const scalarList& nodeAbscissae = abscissae(nodeIndex);

                forAll(node.abscissae(), cmpt)
                {
                    node.abscissae()[cmpt].boundaryFieldRef()[patchi][facei] =
                        nodeAbscissae[cmpt];
                }
            }
        }
    }
}


bool Foam::basicScalarFieldMomentInversion::invertLocalMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    multivariateMomentSet momentsToInvert
    (
        moments.size(),
        momentOrders_,
        moments.supports(),
        smallM0(),
        smallZeta()
    );

    forAll(momentOrders_, momenti)
    {
        const labelList& momentOrder = momentOrders_[momenti];
        momentsToInvert(momentOrder) = moments(momentOrder)[celli];
    }

    if (!momentInverter_().invert(momentsToInvert))
    {
        if (fatalErrorOnFailedRealizabilityTest)
        {
            FatalErrorInFunction
                << "Moment set in cell " << celli << " is not realizable"
                << abort(FatalError);
        }

        return false;
    }

    const mappedScalarList& weights(momentInverter_().weights());
    const mappedList<scalarList>& abscissae(momentInverter_().abscissae());

    forAll(nodes, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        volScalarNode& node = nodes[nodei];

        node.weight()[celli] = weights(nodeIndex);

        const scalarList& nodeAbscissae = abscissae(nodeIndex);

        forAll(node.abscissae(), cmpt)
        {
            node.abscissae()[cmpt][celli] = nodeAbscissae[cmpt];
        }
    }

    return true;
}


void Foam::basicScalarFieldMomentInversion::invert
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    NotImplemented;
}


void Foam::basicScalarFieldMomentInversion::invertBoundaryMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    NotImplemented;
}


bool Foam::basicScalarFieldMomentInversion::invertLocalMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    NotImplemented;

    return false;
}


Foam::scalar Foam::basicScalarFieldMomentInversion::smallM0() const
{
    return momentInverter_().smallM0();
}


Foam::scalar Foam::basicScalarFieldMomentInversion::smallZeta() const
{
    return momentInverter_().smallZeta();
}

// ************************************************************************* //

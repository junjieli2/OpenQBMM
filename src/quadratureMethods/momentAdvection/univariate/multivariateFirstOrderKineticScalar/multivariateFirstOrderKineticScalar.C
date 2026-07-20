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

#include "multivariateFirstOrderKineticScalar.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace univariateAdvection
{
    defineTypeNameAndDebug(multivariateFirstOrderKineticScalar, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentAdvection,
        multivariateFirstOrderKineticScalar,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::univariateAdvection::multivariateFirstOrderKineticScalar::
multivariateFirstOrderKineticScalar
(
    const dictionary& dict,
    const scalarQuadratureApproximation& quadrature,
    const surfaceScalarField& phi,
    const supportType& support
)
:
    univariateMomentAdvection(dict, quadrature, phi, support),
    nodes_(),
    nodesNei_(),
    nodesOwn_(),
    momentsNei_
    (
        name_,
        nMoments_,
        nodesNei_,
        nDimensions_,
        moments_.map(),
        moments_.supports()
    ),
    momentsOwn_
    (
        name_,
        nMoments_,
        nodesOwn_,
        nDimensions_,
        moments_.map(),
        moments_.supports()
    ),
    momentFieldInverter_(),
    nNodes_(quadrature.nodes().size())
{
    const Map<label> nodeMap(quadrature.nodes().map());

    nodes_.reset(new mappedPtrList<volScalarNode>(nNodes_, nodeMap));
    nodesNei_.reset(new mappedPtrList<surfaceScalarNode>(nNodes_, nodeMap));
    nodesOwn_.reset(new mappedPtrList<surfaceScalarNode>(nNodes_, nodeMap));

    mappedPtrList<volScalarNode>& nodes = nodes_();
    mappedPtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    mappedPtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    PtrList<dimensionSet> abscissaDimensions(nDimensions_);
    const labelList zeroOrder(nDimensions_, 0);

    forAll(abscissaDimensions, dimi)
    {
        labelList firstOrder(zeroOrder);
        firstOrder[dimi] = 1;

        abscissaDimensions.set
        (
            dimi,
            new dimensionSet
            (
                moments_(firstOrder).dimensions()/moments_(zeroOrder).dimensions()
            )
        );
    }

    forAll(nodes, nodei)
    {
        const word nodeName
        (
            "nodeAdvection"
          + mappedPtrList<label>::listToWord(quadrature.nodeIndexes()[nodei])
        );

        nodes.set
        (
            nodei,
            new volScalarNode
            (
                nodeName,
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                abscissaDimensions
            )
        );

        nodesNei.set
        (
            nodei,
            new surfaceScalarNode
            (
                nodeName + "Nei",
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                abscissaDimensions
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceScalarNode
            (
                nodeName + "Own",
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                abscissaDimensions
            )
        );
    }

    forAll(momentsNei_, momenti)
    {
        momentsNei_.set
        (
            momenti,
            new surfaceScalarMoment
            (
                name_,
                moments_[momenti].cmptOrders(),
                nodesNei_,
                fvc::interpolate(moments_[momenti]),
                "Nei"
            )
        );

        momentsOwn_.set
        (
            momenti,
            new surfaceScalarMoment
            (
                name_,
                moments_[momenti].cmptOrders(),
                nodesOwn_,
                fvc::interpolate(moments_[momenti]),
                "Own"
            )
        );
    }

    momentFieldInverter_.reset
    (
        new basicScalarFieldMomentInversion
        (
            dict,
            moments_[0].mesh(),
            quadrature.momentOrders(),
            quadrature.nodeIndexes(),
            labelList(1, -1)
        )
    );
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::univariateAdvection::multivariateFirstOrderKineticScalar::
~multivariateFirstOrderKineticScalar()
{}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::univariateAdvection::multivariateFirstOrderKineticScalar::
interpolateNodes()
{
    const mappedPtrList<volScalarNode>& nodes = nodes_();
    mappedPtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    mappedPtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    IStringStream weightOwnLimiter("upwind");
    IStringStream abscissaOwnLimiter("upwind");

    tmp<surfaceInterpolationScheme<scalar>> weightOwnScheme
    (
        fvc::scheme<scalar>(own_, weightOwnLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> abscissaOwnScheme
    (
        fvc::scheme<scalar>(own_, abscissaOwnLimiter)
    );

    IStringStream weightNeiLimiter("upwind");
    IStringStream abscissaNeiLimiter("upwind");

    tmp<surfaceInterpolationScheme<scalar>> weightNeiScheme
    (
        fvc::scheme<scalar>(nei_, weightNeiLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> abscissaNeiScheme
    (
        fvc::scheme<scalar>(nei_, abscissaNeiLimiter)
    );

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        surfaceScalarNode& nodeNei = nodesNei[nodei];
        surfaceScalarNode& nodeOwn = nodesOwn[nodei];

        nodeOwn.weight() = weightOwnScheme().interpolate(node.weight());
        nodeNei.weight() = weightNeiScheme().interpolate(node.weight());

        forAll(node.abscissae(), cmpt)
        {
            nodeOwn.abscissae()[cmpt] =
                abscissaOwnScheme().interpolate(node.abscissae()[cmpt]);

            nodeNei.abscissae()[cmpt] =
                abscissaNeiScheme().interpolate(node.abscissae()[cmpt]);
        }
    }
}


Foam::scalar
Foam::univariateAdvection::multivariateFirstOrderKineticScalar::
realizableCo() const
{
    return 1.0;
}


void Foam::univariateAdvection::multivariateFirstOrderKineticScalar::update()
{
    momentFieldInverter_().invert(moments_, nodes_());
    interpolateNodes();
    momentsNei_.update();
    momentsOwn_.update();

    const dimensionedScalar zeroPhi("zero", phi_.dimensions(), Zero);

    forAll(divMoments_, momenti)
    {
        divMoments_[momenti] =
            fvc::surfaceIntegrate
            (
                momentsNei_[momenti]*min(phi_, zeroPhi)
              + momentsOwn_[momenti]*max(phi_, zeroPhi)
            );
    }
}

// ************************************************************************* //

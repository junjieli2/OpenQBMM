/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenQBMM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenQBMM and is distributed under the GNU General
    Public License, version 3 or later.

Class
    Foam::PDFTransportModels::populationBalanceModels::
    crysDualImpurityPopulationBalance

Description
    Reduced three-coordinate crystal population balance in one
    volume-equivalent size L and two normalized incorporated-impurity masses
    z1 and z2. Each impurity adsorbs on its own Langmuir isotherm; the two
    Kubota--Mullin inhibition terms add and reduce one size-growth rate,
    GL = G0*clip(1 - alpha1*theta1 - alpha2*theta2, 0, 1), with
    alphar = inhibitionScaler*betar/(T*sigma). Crystal shape is not resolved.

\*---------------------------------------------------------------------------*/

#ifndef crysDualImpurityPopulationBalance_H
#define crysDualImpurityPopulationBalance_H

#include "PDFTransportModel.H"
#include "populationBalanceModel.H"
#include "quadratureApproximations.H"
#include "univariateMomentAdvection.H"
#include "realizableOdeSolver.H"
#include "nucleationModel.H"
#include "diffusionModel.H"

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

class crysDualImpurityPopulationBalance
:
    public PDFTransportModel,
    public populationBalanceModel,
    private realizableOdeSolver<volScalarMoment, volScalarNode>
{
    typedef realizableOdeSolver<volScalarMoment, volScalarNode> odeType;

    scalarQuadratureApproximation quadrature_;
    autoPtr<univariateMomentAdvection> momentAdvection_;

    Switch nucleation_;
    Switch growth_;
    Switch speciesCoupled_;

    autoPtr<populationBalanceSubModels::nucleationModel> nucleationModel_;
    autoPtr<populationBalanceSubModels::diffusionModel> diffusionModel_;

    const volScalarField& T_;
    const volScalarField& C_;
    const volScalarField& Ci1_;
    const volScalarField& Ci2_;
    const volScalarField& sigma_;

    dimensionedScalar rhop_;
    dimensionedScalar shapeFactor_;
    dimensionedScalar sizeRef_;
    dimensionedScalar mRef_;
    dimensionedScalar rhoSolvent_;

    dimensionedScalar CgGrowth_;
    scalar gGrowth_;
    scalar growthScale_;

    scalar kads01_;
    scalar kdes01_;
    scalar kads02_;
    scalar kdes02_;
    dimensionedScalar deltaGads1_;
    dimensionedScalar deltaGdes1_;
    dimensionedScalar deltaGads2_;
    dimensionedScalar deltaGdes2_;
    dimensionedScalar gasConstant_;
    dimensionedScalar beta1_;
    dimensionedScalar beta2_;
    scalar inhibitionScale1_;
    scalar inhibitionScale2_;

    dimensionedScalar molarMassKDP_;
    dimensionedScalar molarMassCGM1_;
    dimensionedScalar molarMassCGM2_;
    dimensionedScalar Gmin1_;
    dimensionedScalar Gmin2_;
    dimensionedScalar km01_;
    dimensionedScalar km02_;
    scalar Ke1_;
    scalar Ke2_;
    scalar partitionScale1_;
    scalar partitionScale2_;
    dimensionedScalar z1Nucleation_;
    dimensionedScalar z2Nucleation_;

    volScalarField L10_;
    volScalarField L32_;
    volScalarField meanZ1_;
    volScalarField meanZ2_;
    volScalarField solidCGM1Ppm_;
    volScalarField solidCGM2Ppm_;
    volScalarField totalImpurityPpm_;
    volScalarField crystalPurityPercent_;
    volScalarField inhibitionFactor_;
    volScalarField equivalentGrowth_;
    volScalarField SYact_;
    volScalarField SI1act_;
    volScalarField SI2act_;

    scalarField crystallizationSourceScale_;

    void validateConfiguration() const;
    void updateStatistics();
    void updateKineticDiagnostics();
    void calcSpeciesTransfer
    (
        const scalarField& m300Before,
        const scalarField& m010Before,
        const scalarField& m001Before
    );

    scalar massTransferCoefficient
    (
        const scalar growthRate,
        const scalar km0
    ) const;

    void kineticState
    (
        const label celli,
        scalar& p,
        scalar& GL,
        scalar& q1,
        scalar& q2
    ) const;

    scalar crystallizationSource
    (
        const labelList& momentOrder,
        const label celli
    );

    tmp<fvScalarMatrix> implicitMomentSource(const volScalarMoment& moment);
    void explicitMomentSource();
    bool solveMomentSources() const;
    bool solveMomentOde() const;

    virtual void updateCellMomentSource(const label celli);

    virtual scalar cellMomentSource
    (
        const labelList& momentOrder,
        const label celli,
        const scalarQuadratureApproximation& quadrature,
        const label environment
    );

public:

    TypeName("crysDualImpurityPopulationBalance");

    crysDualImpurityPopulationBalance
    (
        const word& name,
        const dictionary& dict,
        const surfaceScalarField& phi
    );

    virtual ~crysDualImpurityPopulationBalance();

    virtual scalar realizableCo() const;
    virtual scalar CoNum() const;
    virtual void solve();
    virtual bool readIfModified();

    void correctKinetics();

    bool speciesCoupled() const
    {
        return speciesCoupled_;
    }

    const volScalarField::Internal& SYact() const
    {
        return SYact_();
    }

    const volScalarField::Internal& SI1act() const
    {
        return SI1act_();
    }

    const volScalarField::Internal& SI2act() const
    {
        return SI2act_();
    }

    const scalarQuadratureApproximation& quadrature() const
    {
        return quadrature_;
    }
};

} // End namespace populationBalanceModels
} // End namespace PDFTransportModels
} // End namespace Foam

#endif

// ************************************************************************* //

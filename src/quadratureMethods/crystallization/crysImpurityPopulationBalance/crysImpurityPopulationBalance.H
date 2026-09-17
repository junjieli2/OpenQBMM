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
    crysImpurityPopulationBalance

Description
    Bivariate crystal population balance in crystal length L and normalized
    incorporated-impurity mass z.  Crystal growth and nucleation reuse the
    existing crystallization sub-models.  Competitive adsorption reduces
    length growth and drives impurity incorporation.

\*---------------------------------------------------------------------------*/

#ifndef crysImpurityPopulationBalance_H
#define crysImpurityPopulationBalance_H

#include "PDFTransportModel.H"
#include "populationBalanceModel.H"
#include "quadratureApproximations.H"
#include "univariateMomentAdvection.H"
#include "realizableOdeSolver.H"
#include "nucleationModel.H"
#include "growthModel.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "diffusionModel.H"
#include "impurityAdsorptionModel.H"

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

class crysImpurityPopulationBalance
:
    public PDFTransportModel,
    public populationBalanceModel,
    private realizableOdeSolver<volScalarMoment, volScalarNode>
{
    typedef realizableOdeSolver<volScalarMoment, volScalarNode> odeType;

    //- Bivariate (L,z) scalar quadrature
    scalarQuadratureApproximation quadrature_;

    //- Physical-space transport of all mixed moments
    autoPtr<univariateMomentAdvection> momentAdvection_;

    Switch nucleation_;
    Switch growth_;
    Switch aggregation_;
    Switch breakup_;
    Switch speciesCoupled_;

    autoPtr<populationBalanceSubModels::nucleationModel> nucleationModel_;
    autoPtr<populationBalanceSubModels::growthModel> growthModel_;
    autoPtr<populationBalanceSubModels::aggregationKernel> aggregationKernel_;
    autoPtr<populationBalanceSubModels::breakupKernel> breakupKernel_;
    autoPtr<populationBalanceSubModels::diffusionModel> diffusionModel_;
    autoPtr<populationBalanceSubModels::impurityAdsorptionModel>
        impurityAdsorptionModel_;

    dimensionedScalar rhop_;
    dimensionedScalar shapeFactor_;
    dimensionedScalar sizeRef_;
    dimensionedScalar mRef_;
    dimensionedScalar eta_;
    dimensionedScalar rhoi_;
    dimensionedScalar surfaceFactor_;
    dimensionedScalar zNucleation_;

    volScalarField L10_;
    volScalarField L32_;
    volScalarField meanZ_;
    volScalarField solidImpurityFraction_;
    volScalarField SYact_;
    volScalarField SIact_;

    //- Per-cell common limiter for all nucleation and growth sources
    scalarField crystallizationSourceScale_;

    //- Fraction of the available solute a single split step may consume
    scalar speciesConsumptionLimit_;

    //- Host and impurity moments at the start of the current split step
    scalarField m30StepStart_;
    scalarField m01StepStart_;

    //- Solute and impurity available to the current split step
    scalarField soluteStepBudget_;
    scalarField impurityStepBudget_;

    void validateConfiguration() const;
    void subModelsPreUpdate();
    void updateStatistics();
    void calcSpeciesTransfer
    (
        const scalarField& m30Before,
        const scalarField& m01Before
    );

    scalar crystallizationSource
    (
        const labelList& momentOrder,
        const label celli,
        const scalarQuadratureApproximation& quadrature
    );

    scalar pureGrowthSource
    (
        const labelList& momentOrder,
        const label celli,
        const scalarQuadratureApproximation& quadrature
    );

    scalar nucleationSource(const labelList& momentOrder, const label celli) const;

    scalar aggregationSource
    (
        const labelList& momentOrder,
        const label celli,
        const scalarQuadratureApproximation& quadrature,
        const label environment
    ) const;

    scalar breakupSource
    (
        const labelList& momentOrder,
        const label celli,
        const scalarQuadratureApproximation& quadrature,
        const label environment
    ) const;

    tmp<fvScalarMatrix> implicitMomentSource(const volScalarMoment& moment);
    void explicitMomentSource();
    bool solveMomentSources() const;
    bool solveMomentOde() const;

    virtual void updateCellMomentSource(const label celli);

    //- Reject a realizable local step that would overrun the species budget
    virtual bool acceptMomentUpdate(const label celli);

    virtual scalar cellMomentSource
    (
        const labelList& momentOrder,
        const label celli,
        const scalarQuadratureApproximation& quadrature,
        const label environment
    );

public:

    TypeName("crysImpurityPopulationBalance");

    crysImpurityPopulationBalance
    (
        const word& name,
        const dictionary& dict,
        const surfaceScalarField& phi
    );

    virtual ~crysImpurityPopulationBalance();

    virtual scalar realizableCo() const;
    virtual scalar CoNum() const;
    virtual void solve();
    virtual bool readIfModified();

    //- Refresh the algebraic adsorption state after C and Ci are solved
    void correctAdsorption();

    bool speciesCoupled() const
    {
        return speciesCoupled_;
    }

    const volScalarField::Internal& SYact() const
    {
        return SYact_();
    }

    const volScalarField::Internal& SIact() const
    {
        return SIact_();
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

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
\*---------------------------------------------------------------------------*/

#include "crysImpurityPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"
#include "EulerDdtScheme.H"
#include "zeroGradientFvPatchField.H"

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(crysImpurityPopulationBalance, 0);

    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        crysImpurityPopulationBalance,
        dictionary
    );
}
}
}


Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::crysImpurityPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    PDFTransportModel(name, dict, phi.mesh()),
    populationBalanceModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    quadrature_
    (
        name,
        phi.mesh(),
        List<supportType>(2, supportType::RPlus)
    ),
    momentAdvection_
    (
        univariateMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            phi,
            supportType::RPlus
        )
    ),
    nucleation_(dict.lookupOrDefault("nucleation", false)),
    growth_(dict.lookupOrDefault("growth", false)),
    aggregation_(dict.lookupOrDefault("aggregation", false)),
    breakup_(dict.lookupOrDefault("breakup", false)),
    speciesCoupled_(dict.lookupOrDefault("speciesCoupled", false)),
    nucleationModel_(),
    growthModel_(),
    aggregationKernel_(),
    breakupKernel_(),
    diffusionModel_
    (
        populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    impurityAdsorptionModel_
    (
        populationBalanceSubModels::impurityAdsorptionModel::New
        (
            dict.subDict("impurityCoeffs"),
            phi.mesh()
        )
    ),
    rhop_
    (
        dimensionedScalar::getOrDefault
        (
            "rhop",
            dict,
            dimDensity,
            1000.0
        )
    ),
    shapeFactor_
    (
        dimensionedScalar::getOrDefault
        (
            "shapeFactor",
            dict,
            dimless,
            constant::mathematical::pi/6.0
        )
    ),
    sizeRef_
    (
        dimensionedScalar::getOrDefault
        (
            "sizeRef",
            dict,
            dimLength,
            1.0e-6
        )
    ),
    mRef_
    (
        "mRef",
        dimMass,
        rhop_.value()*shapeFactor_.value()*pow3(sizeRef_.value())
    ),
    sourceConsistencyATol_
    (
        readScalar(dict.subDict("odeCoeffs").lookup("ATol"))
    ),
    eta_
    (
        dimensionedScalar::getOrDefault
        (
            "eta",
            dict.subDict("impurityCoeffs"),
            dimless,
            0.0
        )
    ),
    rhoi_
    (
        dimensionedScalar::getOrDefault
        (
            "rhoi",
            dict.subDict("impurityCoeffs"),
            dimDensity,
            0.0
        )
    ),
    surfaceFactor_
    (
        dimensionedScalar::getOrDefault
        (
            "surfaceFactor",
            dict.subDict("impurityCoeffs"),
            dimless,
            1.0
        )
    ),
    zNucleation_
    (
        dimensionedScalar::getOrDefault
        (
            "zNucleation",
            dict.subDict("impurityCoeffs"),
            dimless,
            0.0
        )
    ),
    L10_
    (
        IOobject
        (
            "L10.crysImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    L32_
    (
        IOobject
        (
            "L32.crysImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    meanZ_
    (
        IOobject
        (
            "meanZ.crysImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    solidImpurityFraction_
    (
        IOobject
        (
            "solidImpurityFraction.crysImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    SYact_
    (
        IOobject
        (
            "SYact",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    SIact_
    (
        IOobject
        (
            "SIact",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    crystallizationSourceScale_(phi.mesh().nCells(), 1.0),
    speciesConsumptionLimit_
    (
        dict.lookupOrDefault<scalar>("speciesConsumptionLimit", 0.9)
    ),
    m30StepStart_(phi.mesh().nCells(), 0.0),
    m01StepStart_(phi.mesh().nCells(), 0.0),
    soluteStepBudget_(phi.mesh().nCells(), 0.0),
    impurityStepBudget_(phi.mesh().nCells(), 0.0)
{
    if (speciesConsumptionLimit_ <= 0 || speciesConsumptionLimit_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "speciesConsumptionLimit must lie in (0, 1], not "
            << speciesConsumptionLimit_ << exit(FatalIOError);
    }
    if (nucleation_)
    {
        nucleationModel_ = populationBalanceSubModels::nucleationModel::New
        (
            dict.subDict("nucleationModel"),
            phi.mesh()
        );
    }

    if (growth_)
    {
        dictionary growthDict(dict.subDict("growthModel"));
        const word growthType(growthDict.lookup("growthModel"));

        if
        (
            growthType == "coolCrysGrowthDissolution"
         || growthType == "linearEvaporation"
         || growthType == "nonLinearEvaporation"
        )
        {
            FatalIOErrorInFunction(dict)
                << "crysImpurityPopulationBalance supports positive crystal "
                << "growth only; growth model '" << growthType
                << "' is not supported."
                << exit(FatalIOError);
        }

        if
        (
            growthType == "coolCrysAreaInhibitionGrowth"
         && !growthDict.found("m2Name")
        )
        {
            growthDict.set
            (
                "m2Name",
                IOobject::groupName("moment.20", name)
            );
        }

        growthModel_ = populationBalanceSubModels::growthModel::New
        (
            growthDict,
            phi.mesh()
        );
    }

    if (aggregation_)
    {
        aggregationKernel_ =
            populationBalanceSubModels::aggregationKernel::New
            (
                dict.subDict("aggregationKernel"),
                phi.mesh()
            );
    }

    if (breakup_)
    {
        breakupKernel_ = populationBalanceSubModels::breakupKernel::New
        (
            dict.subDict("breakupKernel"),
            phi.mesh()
        );
    }

    validateConfiguration();
    impurityAdsorptionModel_->preUpdate();
    updateStatistics();

    const scalarField initialM30
    (
        quadrature_.moments()(labelList({3, 0})).primitiveField()
    );
    const scalarField initialM01
    (
        quadrature_.moments()(labelList({0, 1})).primitiveField()
    );
    calcSpeciesTransfer(initialM30, initialM01);

    Info<< "Bivariate crystal coordinates: (L,z), mRef = "
        << mRef_.value() << " kg" << nl
        << "Aggregation maps (L1,z1)+(L2,z2) to "
        << "(cbrt(L1^3+L2^3),z1+z2); breakup preserves z per volume."
        << endl;
}


Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::~crysImpurityPopulationBalance()
{}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::validateConfiguration() const
{
    if (quadrature_.nDimensions() != 2)
    {
        FatalErrorInFunction
            << "crysImpurityPopulationBalance requires exactly two internal "
            << "coordinates (L,z), but found " << quadrature_.nDimensions()
            << abort(FatalError);
    }

    const labelList& nNodes = quadrature_.nNodes();
    const label nL = nNodes[0];
    const label nz = nNodes[1];

    if (nNodes.size() != 2 || nL < 2 || nz < 1)
    {
        FatalErrorInFunction
            << "crysImpurityPopulationBalance requires nNodes = [nL, nz] "
            << "with nL >= 2 and nz >= 1 for the bivariate (L,z) CQMOM "
            << "closure. Found nNodes " << nNodes << abort(FatalError);
    }

    // Build the expected CQMOM moment set for nNodes = [nL, nz]:
    //  - pure L moments:  (i, 0), i = 0 .. 2*nL-1   (2*nL moments)
    //  - mixed moments:   (i, j), i = 0 .. nL-1,
    //                                j = 1 .. 2*nz-1  (nL*(2*nz-1) moments)
    //  Total: nL*(2*nz + 1) moments
    const label nExpectedMoments = nL*(2*nz + 1);
    labelListList expectedMoments(nExpectedMoments);
    {
        label idx = 0;
        for (label i = 0; i < 2*nL; ++i)
        {
            expectedMoments[idx++] = labelList({i, 0});
        }
        for (label j = 1; j <= 2*nz - 1; ++j)
        {
            for (label i = 0; i < nL; ++i)
            {
                expectedMoments[idx++] = labelList({i, j});
            }
        }
    }

    const labelListList& orders = quadrature_.momentOrders();

    if (orders.size() != nExpectedMoments)
    {
        FatalErrorInFunction
            << "For nNodes = [" << nL << ", " << nz << "] the CQMOM closure "
            << "requires " << nExpectedMoments << " moments " << expectedMoments
            << ", but found " << orders.size() << ": " << orders
            << abort(FatalError);
    }

    for (label requiredi = 0; requiredi < nExpectedMoments; ++requiredi)
    {
        bool found = false;

        forAll(orders, orderi)
        {
            found = found
                ||
                (
                    orders[orderi].size() == 2
                 && orders[orderi][0] == expectedMoments[requiredi][0]
                 && orders[orderi][1] == expectedMoments[requiredi][1]
                );
        }

        if (!found)
        {
            FatalErrorInFunction
                << "Missing required moment (" << expectedMoments[requiredi][0]
                << ' ' << expectedMoments[requiredi][1] << ") in " << orders
                << abort(FatalError);
        }
    }

    // Verify node indexes are consistent with nNodes = [nL, nz].
    // The node list must contain every (i, j) with i = 0..nL-1, j = 0..nz-1.
    const labelListList& nodeIndexes = quadrature_.nodeIndexes();
    const label nExpectedNodes = nL*nz;

    if (nodeIndexes.size() != nExpectedNodes)
    {
        FatalErrorInFunction
            << "For nNodes = [" << nL << ", " << nz << "] the node list must "
            << "have " << nExpectedNodes << " entries, but found "
            << nodeIndexes.size() << ": " << nodeIndexes << abort(FatalError);
    }

    for (label i = 0; i < nL; ++i)
    {
        for (label j = 0; j < nz; ++j)
        {
            bool found = false;

            forAll(nodeIndexes, nodei)
            {
                found = found
                    ||
                    (
                        nodeIndexes[nodei].size() == 2
                     && nodeIndexes[nodei][0] == i
                     && nodeIndexes[nodei][1] == j
                    );
            }

            if (!found)
            {
                FatalErrorInFunction
                    << "Missing required node (" << i << ' ' << j
                    << ") in node list " << nodeIndexes << abort(FatalError);
            }
        }
    }

    const dictionary& fieldInversionDict =
        quadrature_.subDict("basicScalarMomentInversion");
    const dictionary& advectionDict = quadrature_.subDict("momentAdvection");
    const dictionary& advectionInversionDict =
        advectionDict.subDict("basicScalarMomentInversion");

    if
    (
        word(quadrature_.lookup("fieldMomentInversion"))
            != "basicScalarFieldMomentInversion"
     || word(advectionDict.lookup("univariateMomentAdvection"))
            != "multivariateFirstOrderKineticScalar"
     || word(fieldInversionDict.lookup("type")) != "conditional"
     || word(advectionInversionDict.lookup("type")) != "conditional"
    )
    {
        FatalErrorInFunction
            << "The bivariate crystal model requires "
            << "basicScalarFieldMomentInversion, "
            << "multivariateFirstOrderKineticScalar and conditional "
            << "moment inversion for both cell and face quadratures."
            << abort(FatalError);
    }

    const wordList fieldSupports
    (
        fieldInversionDict.lookup("supports")
    );
    const wordList advectionSupports
    (
        advectionInversionDict.lookup("supports")
    );

    if
    (
        fieldSupports.size() != 2
     || fieldSupports[0] != "RPlus"
     || fieldSupports[1] != "RPlus"
     || advectionSupports.size() != 2
     || advectionSupports[0] != "RPlus"
     || advectionSupports[1] != "RPlus"
    )
    {
        FatalErrorInFunction
            << "Both bivariate moment inversions require supports "
            << "(RPlus RPlus). Found field supports " << fieldSupports
            << " and advection supports " << advectionSupports
            << abort(FatalError);
    }

    for (label dimi = 0; dimi < 2; ++dimi)
    {
        const word quadratureName("basicQuadrature" + Foam::name(dimi));
        const word fieldType
        (
            fieldInversionDict.subDict(quadratureName).lookup
            (
                "univariateMomentInversion"
            )
        );
        const word advectionType
        (
            advectionInversionDict.subDict(quadratureName).lookup
            (
                "univariateMomentInversion"
            )
        );

        if (fieldType != advectionType)
        {
            FatalErrorInFunction
                << "Cell and face moment inversions must use the same "
                << "univariate closure in coordinate " << dimi << ". Found "
                << fieldType << " and " << advectionType
                << abort(FatalError);
        }
    }

    const mappedPtrList<volScalarNode>& nodes = quadrature_.nodes();

    if
    (
        nodes.empty()
     || nodes[0].sizeIndex() != 0
     || !nodes[0].lengthBased()
     || nodes[0].useVolumeFraction()
    )
    {
        FatalErrorInFunction
            << "Coordinate 0 must be a length-based number-density "
            << "coordinate; volume-fraction weights are not supported."
            << abort(FatalError);
    }

    const volScalarMomentFieldSet& moments = quadrature_.moments();
    const labelList o00({0, 0});
    const labelList o10({1, 0});
    const labelList o01({0, 1});

    if (moments(o00).dimensions() != dimless/dimVolume)
    {
        FatalErrorInFunction
            << "moment.00 must have number-density dimensions "
            << dimless/dimVolume << ", found " << moments(o00).dimensions()
            << abort(FatalError);
    }

    if
    (
        moments(o10).dimensions()/moments(o00).dimensions() != dimLength
     || moments(o01).dimensions()/moments(o00).dimensions() != dimless
    )
    {
        FatalErrorInFunction
            << "The internal-coordinate dimensions must be (length, "
            << "dimensionless) for (L,z)." << abort(FatalError);
    }

    if
    (
        rhop_.dimensions() != dimDensity
     || shapeFactor_.dimensions() != dimless
     || sizeRef_.dimensions() != dimLength
     || mRef_.dimensions() != dimMass
     || eta_.dimensions() != dimless
     || rhoi_.dimensions() != dimDensity
     || surfaceFactor_.dimensions() != dimless
     || zNucleation_.dimensions() != dimless
    )
    {
        FatalErrorInFunction
            << "Expected dimensions: rhop/rhoi " << dimDensity
            << ", sizeRef " << dimLength << ", mRef " << dimMass
            << ", and shapeFactor/eta/surfaceFactor/zNucleation "
            << dimless << ". Found rhop=" << rhop_.dimensions()
            << ", shapeFactor=" << shapeFactor_.dimensions()
            << ", sizeRef=" << sizeRef_.dimensions()
            << ", mRef=" << mRef_.dimensions()
            << ", eta=" << eta_.dimensions()
            << ", rhoi=" << rhoi_.dimensions()
            << ", surfaceFactor=" << surfaceFactor_.dimensions()
            << ", zNucleation=" << zNucleation_.dimensions()
            << abort(FatalError);
    }

    if
    (
        rhop_.value() <= 0
     || shapeFactor_.value() <= 0
     || sizeRef_.value() <= 0
     || mRef_.value() <= 0
     || eta_.value() < 0
     || rhoi_.value() < 0
     || surfaceFactor_.value() < 0
     || zNucleation_.value() < 0
    )
    {
        FatalErrorInFunction
            << "rhop, shapeFactor and sizeRef must be positive; eta, rhoi, "
            << "surfaceFactor and zNucleation must be non-negative."
            << abort(FatalError);
    }

    // If breakup is enabled, verify the configured daughter distribution
    // can evaluate the highest moment order that breakupSource will request.
    // breakupSource maps bivariate moment (i,j) to univariate daughter
    // order (i + 3*j) (volume proportional to L^3).  The maximum is for
    // the highest (i,j) with j > 0: (nL-1) + 3*(2*nz-1).
    if (breakup_)
    {
        const label maxDaughterOrder = (nL - 1) + 3*(2*nz - 1);
        const scalar probeAbscissa = sizeRef_.value();
        const scalar probeValue =
            breakupKernel_->daughterMoment
            (
                maxDaughterOrder,
                probeAbscissa,
                true
            );

        if (probeValue != probeValue || Foam::mag(probeValue) > 1e100)
        {
            FatalErrorInFunction
                << "The configured daughter distribution cannot evaluate "
                << "moment order " << maxDaughterOrder << " required by "
                << "breakupSource for nNodes = [" << nL << ", " << nz
                << "].  Probe at abscissa " << probeAbscissa << " returned "
                << probeValue << ".  Choose a daughter distribution that "
                << "supports this order (e.g. symmetricFragmentation or "
                << "fullFragmentation)." << abort(FatalError);
        }
    }
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::subModelsPreUpdate()
{
    impurityAdsorptionModel_->preUpdate();

    if (aggregation_)
    {
        aggregationKernel_->preUpdate();
    }

    if (breakup_)
    {
        breakupKernel_->preUpdate();
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::implicitMomentSource
(
    const volScalarMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::updateCellMomentSource(const label celli)
{
    crystallizationSourceScale_[celli] = 1.0;

    if (!speciesCoupled_ || (!nucleation_ && !growth_))
    {
        return;
    }

    const labelList o30({3, 0});
    const labelList o01({0, 1});

    // Budget remaining in this split step. The source is evaluated at every
    // stage of the local integrator, so subtracting what the step has already
    // produced keeps the limiter cumulative. Limiting only the instantaneous
    // rate against the start-of-step solute lets a sequence of locally legal
    // sub-steps together remove more solute than is present, which is how the
    // solute field reached negative values.
    const scalar hostGenerated =
        rhop_.value()*shapeFactor_.value()
       *max(quadrature_.moments()(o30)[celli] - m30StepStart_[celli], scalar(0));
    const scalar impurityGenerated =
        mRef_.value()
       *max(quadrature_.moments()(o01)[celli] - m01StepStart_[celli], scalar(0));

    const scalar hostAvailable =
        max(soluteStepBudget_[celli] - hostGenerated, scalar(0));
    const scalar impurityAvailable =
        max(impurityStepBudget_[celli] - impurityGenerated, scalar(0));

    const scalar hostRate =
        rhop_.value()*shapeFactor_.value()
       *max(crystallizationSource(o30, celli, quadrature_), scalar(0));
    const scalar impurityRate =
        mRef_.value()
       *max(crystallizationSource(o01, celli, quadrature_), scalar(0));
    const scalar invDt = 1.0/max
    (
        phi_.mesh().time().deltaTValue(),
        scalar(SMALL)
    );

    // Both caps fall continuously to zero as the remaining budget is used up,
    // so an exhausted cell needs no separate branch. Comparing the rates with
    // SMALL would not be safe here in any case: the impurity rate carries the
    // mRef factor and is legitimately far below SMALL.
    scalar scale = 1.0;

    if (hostRate > 0.0)
    {
        scale = min(scale, hostAvailable*invDt/hostRate);
    }

    if (impurityRate > 0.0)
    {
        scale = min(scale, impurityAvailable*invDt/impurityRate);
    }

    crystallizationSourceScale_[celli] = max(scale, scalar(0));
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::acceptMomentUpdate(const label celli)
{
    if (!speciesCoupled_ || (!nucleation_ && !growth_))
    {
        return true;
    }

    const labelList o30({3, 0});
    const labelList o01({0, 1});

    const scalar hostGenerated =
        rhop_.value()*shapeFactor_.value()
       *max(quadrature_.moments()(o30)[celli] - m30StepStart_[celli], scalar(0));
    const scalar impurityGenerated =
        mRef_.value()
       *max(quadrature_.moments()(o01)[celli] - m01StepStart_[celli], scalar(0));

    if
    (
        hostGenerated > soluteStepBudget_[celli]
     || impurityGenerated > impurityStepBudget_[celli]
    )
    {
        Pout<< "Moment source rejected in cell " << celli
            << ": host=" << hostGenerated
            << "/" << soluteStepBudget_[celli]
            << ", impurity=" << impurityGenerated
            << "/" << impurityStepBudget_[celli] << endl;
    }
    return
        hostGenerated <= soluteStepBudget_[celli]
     && impurityGenerated <= impurityStepBudget_[celli];
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::pureGrowthSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    return growthModel_->phaseSpaceConvection
    (
        momentOrder,
        celli,
        quadrature
    );
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::nucleationSource
(
    const labelList& momentOrder,
    const label celli
) const
{
    if (!nucleation_)
    {
        return 0.0;
    }

    const scalar zFactor = momentOrder[1] == 0
        ? 1.0
        : pow(zNucleation_.value(), momentOrder[1]);

    return nucleationModel_->nucleationSource(momentOrder[0], celli)*zFactor;
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::aggregationSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
) const
{
    scalar source = 0.0;
    const mappedPtrList<volScalarNode>& nodes = quadrature.nodes();
    const label i = momentOrder[0];
    const label j = momentOrder[1];

    // These two moments are invariants of the selected coalescence map.
    if ((i == 3 && j == 0) || (i == 0 && j == 1))
    {
        return 0.0;
    }

    forAll(nodes, node1i)
    {
        const volScalarNode& node1 = nodes[node1i];
        const scalar L1 = max(node1.abscissae()[0][celli], scalar(0));
        const scalar z1 = max(node1.abscissae()[1][celli], scalar(0));
        const scalar n1 = node1.numberDensity
        (
            celli,
            node1.weight()[celli],
            L1
        );
        const scalar d1 = node1.d(celli, L1);

        forAll(nodes, node2i)
        {
            const volScalarNode& node2 = nodes[node2i];
            const scalar L2 = max(node2.abscissae()[0][celli], scalar(0));
            const scalar z2 = max(node2.abscissae()[1][celli], scalar(0));
            const scalar n2 = node2.numberDensity
            (
                celli,
                node2.weight()[celli],
                L2
            );
            const scalar d2 = node2.d(celli, L2);
            const scalar Lnew = pow(pow3(L1) + pow3(L2), 1.0/3.0);
            const scalar znew = z1 + z2;

            const scalar birth = pow(Lnew, i)*pow(znew, j);
            const scalar loss1 = pow(L1, i)*pow(z1, j);
            const scalar loss2 = pow(L2, i)*pow(z2, j);

            source +=
                0.5*n1*n2
               *aggregationKernel_->Ka(d1, d2, Zero, celli, environment)
               *(birth - loss1 - loss2);
        }
    }

    return source;
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::breakupSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
) const
{
    scalar source = 0.0;
    const mappedPtrList<volScalarNode>& nodes = quadrature.nodes();
    const label i = momentOrder[0];
    const label j = momentOrder[1];
    const label daughterOrder = i + 3*j;

    // Volume and incorporated-impurity mass are exactly conserved.
    if ((i == 3 && j == 0) || (i == 0 && j == 1))
    {
        return 0.0;
    }

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        const scalar L = max(node.abscissae()[0][celli], scalar(0));
        const scalar z = max(node.abscissae()[1][celli], scalar(0));
        const scalar n = node.numberDensity
        (
            celli,
            node.weight()[celli],
            L
        );

        scalar daughterTerm = breakupKernel_->daughterMoment
        (
            daughterOrder,
            L,
            true
        );

        if (j > 0)
        {
            daughterTerm /= pow(max(L, scalar(SMALL)), 3*j);
        }

        source +=
            n*breakupKernel_->Kb(L, celli, environment)*pow(z, j)
           *(daughterTerm - pow(L, i));
    }

    return source;
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::crystallizationSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar source = nucleationSource(momentOrder, celli);

    if (growth_)
    {
        const label i = momentOrder[0];
        const label j = momentOrder[1];
        const scalar theta = min
        (
            max(impurityAdsorptionModel_->theta(celli), scalar(0)),
            scalar(1)
        );
        const scalar cleanFraction = 1.0 - theta;

        if (i > 0)
        {
            source += cleanFraction*max
            (
                pureGrowthSource(momentOrder, celli, quadrature),
                scalar(0)
            );
        }

        if (j > 0 && eta_.value() > 0 && rhoi_.value() > 0)
        {
            labelList shiftedOrder(momentOrder);
            shiftedOrder[0] += 3;
            shiftedOrder[1] -= 1;

            const scalar incorporationCoeff =
                eta_.value()*rhoi_.value()*surfaceFactor_.value()
               /mRef_.value();

            source +=
                j*incorporationCoeff*theta*cleanFraction
               *max
                (
                    pureGrowthSource(shiftedOrder, celli, quadrature),
                    scalar(0)
                )
               /scalar(i + 3);
        }
    }

    return source;
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source =
        crystallizationSourceScale_[celli]
       *crystallizationSource(momentOrder, celli, quadrature);

    if (aggregation_)
    {
        source += aggregationSource
        (
            momentOrder,
            celli,
            quadrature,
            environment
        );
    }

    if (breakup_)
    {
        source += breakupSource
        (
            momentOrder,
            celli,
            quadrature,
            environment
        );
    }

    return source;
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::solveMomentSources() const
{
    return
        (nucleation_ || growth_ || aggregation_ || breakup_)
     && odeType::solveSources_;
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::solveMomentOde() const
{
    return odeType::solveOde_;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::updateStatistics()
{
    const volScalarMomentFieldSet& moments = quadrature_.moments();
    const volScalarMoment& m00 = moments(labelList({0, 0}));
    const volScalarMoment& m10 = moments(labelList({1, 0}));
    const volScalarMoment& m20 = moments(labelList({2, 0}));
    const volScalarMoment& m30 = moments(labelList({3, 0}));
    const volScalarMoment& m01 = moments(labelList({0, 1}));

    forAll(L10_, celli)
    {
        L10_[celli] = m10[celli]/max(m00[celli], scalar(SMALL));
        L32_[celli] = m30[celli]/max(m20[celli], scalar(SMALL));
        meanZ_[celli] = m01[celli]/max(m00[celli], scalar(SMALL));

        const scalar crystalMass =
            rhop_.value()*shapeFactor_.value()*max(m30[celli], scalar(0));
        const scalar impurityMass = mRef_.value()*max(m01[celli], scalar(0));

        solidImpurityFraction_[celli] =
            impurityMass/max(crystalMass + impurityMass, scalar(SMALL));
    }

    L10_.correctBoundaryConditions();
    L32_.correctBoundaryConditions();
    meanZ_.correctBoundaryConditions();
    solidImpurityFraction_.correctBoundaryConditions();
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::calcSpeciesTransfer
(
    const scalarField& m30Before,
    const scalarField& m01Before
)
{
    SYact_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);
    SIact_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);

    if
    (
        !speciesCoupled_
     || !solveMomentSources()
     || (!nucleation_ && !growth_)
    )
    {
        SYact_.correctBoundaryConditions();
        SIact_.correctBoundaryConditions();
        return;
    }

    const fvMesh& mesh = phi_.mesh();
    const scalar invDt = 1.0/max(mesh.time().deltaTValue(), scalar(SMALL));
    const volScalarMoment& m30 =
        quadrature_.moments()(labelList({3, 0}));
    const volScalarMoment& m01 =
        quadrature_.moments()(labelList({0, 1}));
    const scalar maximumImpurityHostRatio =
        eta_.value()*rhoi_.value()*surfaceFactor_.value()
       /(3.0*rhop_.value()*shapeFactor_.value());

    forAll(SYact_, celli)
    {
        const scalar deltaM30 = m30[celli] - m30Before[celli];
        const scalar deltaM01 = m01[celli] - m01Before[celli];
        const scalar hostMassDelta =
            rhop_.value()*shapeFactor_.value()*deltaM30;
        const scalar impurityMassDelta = mRef_.value()*deltaM01;
        const scalar hostMassTolerance =
            100.0*SMALL*max(soluteStepBudget_[celli], scalar(VSMALL));
        const scalar impurityMassTolerance =
            100.0*SMALL*max(impurityStepBudget_[celli], scalar(VSMALL));

        const scalar hostInventory =
            rhop_.value()*shapeFactor_.value()
           *max(m30Before[celli], scalar(0));
        const scalar impurityInventory =
            mRef_.value()*max(m01Before[celli], scalar(0));
        const scalar physicalScale = max
        (
            max
            (
                impurityInventory,
                maximumImpurityHostRatio*hostInventory
            ),
            max
            (
                max(impurityMassDelta, scalar(0)),
                maximumImpurityHostRatio*max(hostMassDelta, scalar(0))
            )
        );
        const scalar consistencyTolerance = max
        (
            10.0*mRef_.value()*sourceConsistencyATol_,
            1.0e-6*max(physicalScale, scalar(VSMALL))
        );

        if
        (
            hostMassDelta < -hostMassTolerance
         || impurityMassDelta < -impurityMassTolerance
        )
        {
            FatalErrorInFunction
                << "Negative crystallization source increment in cell " << celli
                << ": host mass=" << hostMassDelta
                << ", impurity mass=" << impurityMassDelta
                << abort(FatalError);
        }

        if
        (
            impurityMassDelta
          > maximumImpurityHostRatio*max(hostMassDelta, scalar(0))
          + consistencyTolerance
        )
        {
            FatalErrorInFunction
                << "Impurity incorporation exceeds the analytic source "
                << "ceiling in cell " << celli
                << ": host mass=" << hostMassDelta
                << ", impurity mass=" << impurityMassDelta
                << ", maximum ratio=" << maximumImpurityHostRatio
                << ", tolerance=" << consistencyTolerance
                << abort(FatalError);
        }
        SYact_[celli] =
            rhop_.value()*shapeFactor_.value()
           *max
            (
                deltaM30 > 0.0 ? deltaM30*invDt : 0.0,
                scalar(0)
            );
        SIact_[celli] =
            mRef_.value()
           *max
            (
                deltaM01 > 0.0 ? deltaM01*invDt : 0.0,
                scalar(0)
            );
    }

    SYact_.correctBoundaryConditions();
    SIact_.correctBoundaryConditions();

    // How close the split step came to the positivity limit. A value near one
    // means the step size, not the kinetics, is setting the transfer and the
    // result is no longer a converged solution of the model.
    const scalar dt = mesh.time().deltaTValue();
    scalar maxHostFraction = 0.0;
    scalar maxImpurityFraction = 0.0;

    forAll(SYact_, celli)
    {
        const scalar hostAvailable =
            max(impurityAdsorptionModel_->solute()[celli], scalar(0));
        const scalar impurityAvailable =
            max(impurityAdsorptionModel_->impurity()[celli], scalar(0));

        if (hostAvailable > SMALL)
        {
            maxHostFraction =
                max(maxHostFraction, SYact_[celli]*dt/hostAvailable);
        }

        if (impurityAvailable > SMALL)
        {
            maxImpurityFraction =
                max(maxImpurityFraction, SIact_[celli]*dt/impurityAvailable);
        }
    }

    reduce(maxHostFraction, maxOp<scalar>());
    reduce(maxImpurityFraction, maxOp<scalar>());

    Info<< "Species step consumption fraction (limit "
        << speciesConsumptionLimit_ << ") host = " << maxHostFraction
        << ", impurity = " << maxImpurityFraction << endl;
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::realizableCo() const
{
    return momentAdvection_->realizableCo();
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::CoNum() const
{
    return 0.0;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::solve()
{
    // This auxiliary field is not a registered geometric field, so keep it
    // synchronized explicitly if a dynamic mesh changes the cell count.
    crystallizationSourceScale_.setSize(phi_.mesh().nCells(), 1.0);
    crystallizationSourceScale_ = 1.0;
    m30StepStart_.setSize(phi_.mesh().nCells(), 0.0);
    m01StepStart_.setSize(phi_.mesh().nCells(), 0.0);
    soluteStepBudget_.setSize(phi_.mesh().nCells(), 0.0);
    impurityStepBudget_.setSize(phi_.mesh().nCells(), 0.0);

    subModelsPreUpdate();
    momentAdvection_->update();

    // Source integration is operator-split from physical-space transport.
    // Multi-step ddt schemes would re-use the already applied source jump on
    // the following time level, so advance these moments with Euler as well.
    fv::EulerDdtScheme<scalar> momentDdt(phi_.mesh());

    forAll(quadrature_.moments(), momenti)
    {
        volScalarMoment& moment = quadrature_.moments()[momenti];

        fvScalarMatrix momentEqn
        (
            momentDdt.fvmDdt(moment)
          + momentAdvection_->divMoments()[momenti]
         == implicitMomentSource(moment)
        );

        if (diffusionModel_->type() != "none")
        {
            momentEqn.relax();
        }

        momentEqn.solve();
        moment.correctBoundaryConditions();
    }

    // Keep the integrated moments as the primary state. Projecting the
    // quadrature back onto the moments after transport changes the host
    // moments but not the impurity moments, which inflates the accumulated
    // incorporation ratio beyond the analytic growth ceiling. Only cells whose
    // inversion fails are recovered from the last valid nodes.
    const label nRecovered = quadrature_.updateQuadraturePreservingMoments();

    if (nRecovered > 0)
    {
        Info<< "Recovered " << nRecovered
            << " cells from failed moment inversion after transport" << nl;
    }

    const scalarField m30Before
    (
        quadrature_.moments()(labelList({3, 0})).primitiveField()
    );
    const scalarField m01Before
    (
        quadrature_.moments()(labelList({0, 1})).primitiveField()
    );

    // Freeze the budgets for the source sub-step. The solute and impurity
    // fields are only advanced by their own equations after this point, so
    // these are the amounts actually available to the reaction step.
    m30StepStart_ = m30Before;
    m01StepStart_ = m01Before;

    if (speciesCoupled_)
    {
        const volScalarField& solute = impurityAdsorptionModel_->solute();
        const volScalarField& impurity = impurityAdsorptionModel_->impurity();

        forAll(soluteStepBudget_, celli)
        {
            soluteStepBudget_[celli] =
                speciesConsumptionLimit_*max(solute[celli], scalar(0));
            impurityStepBudget_[celli] =
                speciesConsumptionLimit_*max(impurity[celli], scalar(0));
        }
    }

    if (solveMomentSources())
    {
        explicitMomentSource();
    }

    updateStatistics();
    calcSpeciesTransfer(m30Before, m01Before);
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_.subDict(type() + "Coeffs")
    );
    return true;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysImpurityPopulationBalance::correctAdsorption()
{
    impurityAdsorptionModel_->correct();
}

// ************************************************************************* //

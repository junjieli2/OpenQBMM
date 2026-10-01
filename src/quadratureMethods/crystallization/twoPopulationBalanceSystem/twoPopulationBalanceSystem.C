#include "twoPopulationBalanceSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "EulerDdtScheme.H"
#include "zeroGradientFvPatchField.H"

#include <algorithm>
#include <cmath>

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(twoPopulationBalanceSystem, 0);

    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        twoPopulationBalanceSystem,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace
{
Foam::Switch observedMomentEnabled(const Foam::dictionary& dict)
{
    return
        dict.found("observationModel")
      ? dict.subDict("observationModel").lookupOrDefault<Foam::Switch>
        (
            "enabled",
            true
        )
      : Foam::Switch(false);
}


Foam::dimensionedScalar detectionL50(const Foam::dictionary& dict)
{
    const Foam::dimensionedScalar defaultValue
    (
        "L50",
        Foam::dimLength,
        50.0e-6
    );

    return
        dict.found("observationModel")
      ? dict.subDict("observationModel").lookupOrDefault
        (
            "L50",
            defaultValue
        )
      : defaultValue;
}


Foam::scalar detectionLogWidth(const Foam::dictionary& dict)
{
    return
        dict.found("observationModel")
      ? dict.subDict("observationModel").lookupOrDefault<Foam::scalar>
        (
            "logWidth",
            0.25
        )
      : 0.25;
}
}


Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::twoPopulationBalanceSystem
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    populationBalanceModel(name, dict, phi),
    singleCrystalDict_(dict.subDict("singleCrystal")),
    agglomerateDict_(dict.subDict("agglomerate")),
    ssDict_(dict.subDict("aggregation").subDict("ss")),
    saDict_(dict.subDict("aggregation").subDict("sa")),
    aaDict_(dict.subDict("aggregation").subDict("aa")),
    singleCrystalQuadrature_
    (
        "singleCrystal",
        phi.mesh(),
        List<supportType>(1, supportType::RPlus)
    ),
    agglomerateQuadrature_
    (
        "agglomerate",
        phi.mesh(),
        List<supportType>(1, supportType::RPlus)
    ),
    singleCrystalAdvection_
    (
        univariateMomentAdvection::New
        (
            singleCrystalQuadrature_.subDict("momentAdvection"),
            singleCrystalQuadrature_,
            phi,
            supportType::RPlus
        )
    ),
    agglomerateAdvection_
    (
        univariateMomentAdvection::New
        (
            agglomerateQuadrature_.subDict("momentAdvection"),
            agglomerateQuadrature_,
            phi,
            supportType::RPlus
        )
    ),
    nucleation_(singleCrystalDict_.lookupOrDefault("nucleation", false)),
    singleGrowth_(singleCrystalDict_.lookupOrDefault("growth", false)),
    agglomerateGrowth_(agglomerateDict_.lookupOrDefault("growth", false)),
    ssEnabled_(ssDict_.lookupOrDefault("enabled", false)),
    saEnabled_(saDict_.lookupOrDefault("enabled", false)),
    aaEnabled_(aaDict_.lookupOrDefault("enabled", false)),
    speciesCoupled_(dict.lookupOrDefault("speciesCoupled", false)),
    agglomerateBreakup_(agglomerateDict_.lookupOrDefault("breakup", false)),
    singleCrystalDaughterFraction_
    (
        agglomerateDict_.lookupOrDefault<scalar>
        (
            "singleCrystalDaughterFraction", 0.0
        )
    ),
    nucleationModel_(),
    singleGrowthModel_(),
    agglomerateGrowthModel_(),
    singleDiffusionModel_
    (
        populationBalanceSubModels::diffusionModel::New
        (
            singleCrystalDict_.subDict("diffusionModel")
        )
    ),
    agglomerateDiffusionModel_
    (
        populationBalanceSubModels::diffusionModel::New
        (
            agglomerateDict_.subDict("diffusionModel")
        )
    ),
    crossGrowthRate_
    (
        IOobject
        (
            "growthRate.sa",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    ssKernel_(),
    saKernel_(),
    aaKernel_(),
    agglomerateBreakupKernel_(),
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
    ATol_(readScalar(dict.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(dict.subDict("odeCoeffs").lookup("RTol"))),
    fac_(readScalar(dict.subDict("odeCoeffs").lookup("fac"))),
    facMin_(readScalar(dict.subDict("odeCoeffs").lookup("facMin"))),
    facMax_(readScalar(dict.subDict("odeCoeffs").lookup("facMax"))),
    minLocalDt_
    (
        readScalar(dict.subDict("odeCoeffs").lookup("minLocalDt"))
    ),
    solveSources_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("solveSources", true)
    ),
    solveOde_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("solveOde", true)
    ),
    localDt_
    (
        IOobject
        (
            "twoPopulationOde.localDt",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        phi.mesh().time().deltaT()
    ),
    localDtAdjustments_(0),
    totalM2_
    (
        IOobject
        (
            "moment.2.twoPopulation",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar
        (
            "zero",
            singleCrystalQuadrature_.moments()
            (
                labelList(1, 2)
            ).dimensions(),
            0.0
        ),
        fvPatchFieldBase::zeroGradientType()
    ),
    singleL10_
    (
        IOobject
        (
            "L10.singleCrystal",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    singleL32_
    (
        IOobject
        (
            "L32.singleCrystal",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    agglomerateL10_
    (
        IOobject
        (
            "L10.agglomerate",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    agglomerateL32_
    (
        IOobject
        (
            "L32.agglomerate",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    totalL10_
    (
        IOobject
        (
            "L10.total",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    totalSolidMass_
    (
        IOobject
        (
            "solidMassDensity.twoPopulation",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    SYact_
    (
        IOobject
        (
            "SYact.twoPopulation",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    observedMoment_(observedMomentEnabled(dict)),
    detectionL50_(detectionL50(dict)),
    detectionLogWidth_(detectionLogWidth(dict)),
    observedMomentOrders_(),
    observedMoments_(),
    observedD32_()
{
    updateTotalM2();

    if (singleGrowth_)
    {
        dictionary& growthDict = singleCrystalDict_.subDict("growthModel");
        growthDict.set("growthRateField", word("growthRate.singleCrystal"));
        singleGrowthModel_ = populationBalanceSubModels::growthModel::New
        (
            growthDict,
            phi.mesh()
        );
    }

    if (agglomerateGrowth_)
    {
        dictionary& growthDict = agglomerateDict_.subDict("growthModel");
        growthDict.set("growthRateField", word("growthRate.agglomerate"));
        agglomerateGrowthModel_ =
            populationBalanceSubModels::growthModel::New
            (
                growthDict,
                phi.mesh()
            );
    }

    if (nucleation_)
    {
        dictionary& nucleationDict =
            singleCrystalDict_.subDict("nucleationModel");

        const word nucleationType
        (
            nucleationDict.lookup("nucleationModel")
        );

        if
        (
            nucleationType == "areaNucleation"
         && !nucleationDict.found("m2Name")
        )
        {
            nucleationDict.set("m2Name", totalM2_.name());
        }
        else if (!nucleationDict.found("momentGroup"))
        {
            nucleationDict.set("momentGroup", word("singleCrystal"));
        }

        nucleationModel_ = populationBalanceSubModels::nucleationModel::New
        (
            nucleationDict,
            phi.mesh()
        );
    }

    if (ssEnabled_)
    {
        if (!ssDict_.found("growthRateField"))
        {
            ssDict_.set
            (
                "growthRateField",
                word("growthRate.singleCrystal")
            );
        }
        if (!ssDict_.found("aggregationRateField"))
        {
            ssDict_.set("aggregationRateField", word("aggregationRate.ss"));
        }
        if (!ssDict_.found("L10Field"))
        {
            ssDict_.set("L10Field", singleL10_.name());
        }
        ssKernel_ = populationBalanceSubModels::aggregationKernel::New
        (
            ssDict_,
            phi.mesh()
        );
    }

    if (saEnabled_)
    {
        saDict_.set("growthRateField", word("growthRate.sa"));
        if (!saDict_.found("aggregationRateField"))
        {
            saDict_.set("aggregationRateField", word("aggregationRate.sa"));
        }
        if (!saDict_.found("L10Field"))
        {
            saDict_.set("L10Field", totalL10_.name());
        }
        saKernel_ = populationBalanceSubModels::aggregationKernel::New
        (
            saDict_,
            phi.mesh()
        );
    }

    if (aaEnabled_)
    {
        if (!aaDict_.found("growthRateField"))
        {
            aaDict_.set
            (
                "growthRateField",
                word("growthRate.agglomerate")
            );
        }
        if (!aaDict_.found("aggregationRateField"))
        {
            aaDict_.set("aggregationRateField", word("aggregationRate.aa"));
        }
        if (!aaDict_.found("L10Field"))
        {
            aaDict_.set("L10Field", agglomerateL10_.name());
        }
        aaKernel_ = populationBalanceSubModels::aggregationKernel::New
        (
            aaDict_,
            phi.mesh()
        );
    }

    if (agglomerateBreakup_)
    {
        agglomerateBreakupKernel_ = populationBalanceSubModels::breakupKernel::New
        (
            agglomerateDict_.subDict("breakupKernel"), phi.mesh()
        );
    }

    validateConfiguration();

    if (observedMoment_)
    {
        if (detectionL50_.value() <= 0.0)
        {
            FatalIOErrorInFunction(dict)
                << "observationModel.L50 must be positive, found "
                << detectionL50_ << exit(FatalIOError);
        }

        if (detectionLogWidth_ <= SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "observationModel.logWidth must be positive, found "
                << detectionLogWidth_ << exit(FatalIOError);
        }

        const scalarQuadratureApproximation* quadratures[2] =
        {
            &singleCrystalQuadrature_,
            &agglomerateQuadrature_
        };

        for (label populationi = 0; populationi < 2; ++populationi)
        {
            const labelListList& momentOrders =
                quadratures[populationi]->momentOrders();

            forAll(momentOrders, momenti)
            {
                const label order = momentOrders[momenti][0];

                if (!observedMomentOrders_.found(order))
                {
                    observedMomentOrders_.append(order);
                }
            }
        }

        std::sort
        (
            observedMomentOrders_.begin(),
            observedMomentOrders_.end()
        );

        observedMoments_.setSize(observedMomentOrders_.size());

        forAll(observedMomentOrders_, momenti)
        {
            const label order = observedMomentOrders_[momenti];
            const label singleIndex =
                momentIndex(singleCrystalQuadrature_, order);
            const label agglomerateIndex =
                momentIndex(agglomerateQuadrature_, order);

            const dimensionSet& dimensions =
                singleIndex >= 0
              ? singleCrystalQuadrature_.moments()[singleIndex].dimensions()
              : agglomerateQuadrature_.moments()[agglomerateIndex].dimensions();

            observedMoments_.set
            (
                momenti,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            word("observedMoment.") + Foam::name(order),
                            name
                        ),
                        phi.mesh().time().timeName(),
                        phi.mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    phi.mesh(),
                    dimensionedScalar("zero", dimensions, 0.0),
                    fvPatchFieldBase::zeroGradientType()
                )
            );
        }

        observedD32_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("observedD32", name),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimLength, 0.0),
                fvPatchFieldBase::zeroGradientType()
            )
        );

        Info<< "Camera observation model enabled for the merged populations:"
            << nl
            << "  L50 = " << detectionL50_.value() << " m" << nl
            << "  logWidth = " << detectionLogWidth_ << endl;
    }

    preUpdate();
    updateStatistics();
    updateObservedMoments();

    scalarField totalM3Before
    (
        singleCrystalQuadrature_.moments()
        (
            labelList(1, 3)
        ).primitiveField()
      + agglomerateQuadrature_.moments()
        (
            labelList(1, 3)
        ).primitiveField()
    );
    calculateSpeciesTransfer(totalM3Before);

    Info<< "Two-population crystal PBE: nucleation enters singleCrystal; "
        << "ss, sa and aa collision products enter agglomerate. "
        << "Coalescence uses Lnew=cbrt(L1^3+L2^3)." << endl;
}


Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::~twoPopulationBalanceSystem()
{}


Foam::label Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::momentIndex
(
    const scalarQuadratureApproximation& quadrature,
    const label order
) const
{
    const labelListList& orders = quadrature.momentOrders();

    forAll(orders, momenti)
    {
        if (orders[momenti].size() == 1 && orders[momenti][0] == order)
        {
            return momenti;
        }
    }

    return -1;
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::validateConfiguration() const
{
    const scalarQuadratureApproximation* quadratures[2] =
    {
        &singleCrystalQuadrature_,
        &agglomerateQuadrature_
    };

    const word names[2] = {"singleCrystal", "agglomerate"};

    for (label populationi = 0; populationi < 2; ++populationi)
    {
        const scalarQuadratureApproximation& quadrature =
            *quadratures[populationi];

        if (quadrature.nDimensions() != 1)
        {
            FatalErrorInFunction
                << names[populationi] << " must be univariate."
                << abort(FatalError);
        }

        for (label order = 0; order <= 3; ++order)
        {
            if (momentIndex(quadrature, order) < 0)
            {
                FatalErrorInFunction
                    << names[populationi]
                    << " quadrature requires moments 0 through 3."
                    << abort(FatalError);
            }
        }

        if
        (
            quadrature.nodes().empty()
         || !quadrature.nodes()[0].lengthBased()
         || quadrature.nodes()[0].sizeIndex() != 0
         || quadrature.nodes()[0].useVolumeFraction()
        )
        {
            FatalErrorInFunction
                << names[populationi]
                << " must use a length-based number-density coordinate."
                << abort(FatalError);
        }

        const volScalarMomentFieldSet& moments = quadrature.moments();
        if
        (
            moments(labelList(1, 0)).dimensions() != dimless/dimVolume
         || moments(labelList(1, 1)).dimensions()
               /moments(labelList(1, 0)).dimensions() != dimLength
        )
        {
            FatalErrorInFunction
                << names[populationi]
                << " moments must be number-density moments of length."
                << abort(FatalError);
        }
    }

    if
    (
        !std::isfinite(singleCrystalDaughterFraction_)
     || singleCrystalDaughterFraction_ < 0
     || singleCrystalDaughterFraction_ > 1
    )
    {
        FatalIOErrorInFunction(agglomerateDict_)
            << "singleCrystalDaughterFraction must be in [0,1]."
            << exit(FatalIOError);
    }

    if (rhop_.value() <= 0 || shapeFactor_.value() <= 0)
    {
        FatalErrorInFunction
            << "rhop and shapeFactor must be positive."
            << abort(FatalError);
    }

    if
    (
        ATol_ <= 0
     || RTol_ <= 0
     || fac_ <= 0
     || facMin_ <= 0
     || facMax_ < facMin_
     || minLocalDt_ <= 0
    )
    {
        FatalErrorInFunction
            << "Invalid odeCoeffs for twoPopulationBalanceSystem."
            << abort(FatalError);
    }

    if
    (
        speciesCoupled_
     && !phi_.mesh().foundObject<volScalarField>("C")
    )
    {
        FatalErrorInFunction
            << "speciesCoupled requires the solute field C."
            << abort(FatalError);
    }
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::updateTotalM2()
{
    const volScalarMoment& singleM2 =
        singleCrystalQuadrature_.moments()(labelList(1, 2));
    const volScalarMoment& agglomerateM2 =
        agglomerateQuadrature_.moments()(labelList(1, 2));

    forAll(totalM2_, celli)
    {
        totalM2_[celli] =
            max(singleM2[celli] + agglomerateM2[celli], scalar(0));
    }

    totalM2_.correctBoundaryConditions();
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::updateCrossGrowthRate()
{
    const fvMesh& mesh = phi_.mesh();
    const volScalarField* singleGrowthField = nullptr;
    const volScalarField* agglomerateGrowthField = nullptr;

    if (mesh.foundObject<volScalarField>("growthRate.singleCrystal"))
    {
        singleGrowthField =
            &mesh.lookupObject<volScalarField>("growthRate.singleCrystal");
    }

    if (mesh.foundObject<volScalarField>("growthRate.agglomerate"))
    {
        agglomerateGrowthField =
            &mesh.lookupObject<volScalarField>("growthRate.agglomerate");
    }

    forAll(crossGrowthRate_, celli)
    {
        const scalar Gs =
            singleGrowthField ? (*singleGrowthField)[celli] : scalar(0);
        const scalar Ga =
            agglomerateGrowthField
          ? (*agglomerateGrowthField)[celli]
          : scalar(0);

        crossGrowthRate_[celli] = 0.5*(Gs + Ga);
    }

    crossGrowthRate_.correctBoundaryConditions();
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::preUpdate()
{
    const labelList firstOrder(1, 1);

    forAll(singleCrystalQuadrature_.moments()[0], celli)
    {
        if (singleGrowth_)
        {
            singleGrowthModel_->phaseSpaceConvection
            (
                firstOrder,
                celli,
                singleCrystalQuadrature_
            );
        }

        if (agglomerateGrowth_)
        {
            agglomerateGrowthModel_->phaseSpaceConvection
            (
                firstOrder,
                celli,
                agglomerateQuadrature_
            );
        }
    }

    updateCrossGrowthRate();

    if (agglomerateBreakup_)
    {
        agglomerateBreakupKernel_->preUpdate();
    }
    if (ssEnabled_)
    {
        ssKernel_->preUpdate();
    }
    if (saEnabled_)
    {
        saKernel_->preUpdate();
    }
    if (aaEnabled_)
    {
        aaKernel_->preUpdate();
    }
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::aggregationKernelValue
(
    const populationBalanceSubModels::aggregationKernel& kernel,
    const volScalarNode& first,
    const volScalarNode& second,
    const label celli
) const
{
    const scalar L1 = max(first.abscissae()[0][celli], scalar(0));
    const scalar L2 = max(second.abscissae()[0][celli], scalar(0));
    const scalar d1 = first.d(celli, L1);
    const scalar d2 = second.d(celli, L2);

    const scalar value = kernel.Ka(d1, d2, Zero, celli, 0);
    return std::isfinite(value) && value > 0 ? value : scalar(0);
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::calculateSources
(
    const label celli,
    scalarList& singleSources,
    scalarList& agglomerateSources
)
{
    const label nSingle = singleCrystalQuadrature_.nMoments();
    const label nAgglomerate = agglomerateQuadrature_.nMoments();
    singleSources.setSize(nSingle, 0.0);
    agglomerateSources.setSize(nAgglomerate, 0.0);

    scalarList singleCrystallization(nSingle, 0.0);
    scalarList agglomerateCrystallization(nAgglomerate, 0.0);
    scalarList singleAggregation(nSingle, 0.0);

    const labelListList& singleOrders =
        singleCrystalQuadrature_.momentOrders();
    const labelListList& agglomerateOrders =
        agglomerateQuadrature_.momentOrders();

    const label singleM2 = momentIndex(singleCrystalQuadrature_, 2);
    const label agglomerateM2 = momentIndex(agglomerateQuadrature_, 2);
    totalM2_[celli] = max
    (
        singleCrystalQuadrature_.moments()[singleM2][celli]
      + agglomerateQuadrature_.moments()[agglomerateM2][celli],
        scalar(0)
    );

    forAll(singleCrystallization, momenti)
    {
        const label order = singleOrders[momenti][0];

        if (nucleation_)
        {
            singleCrystallization[momenti] +=
                nucleationModel_->nucleationSource(order, celli);
        }

        if (singleGrowth_)
        {
            singleCrystallization[momenti] +=
                singleGrowthModel_->phaseSpaceConvection
                (
                    singleOrders[momenti],
                    celli,
                    singleCrystalQuadrature_
                );
        }
    }

    forAll(agglomerateCrystallization, momenti)
    {
        if (agglomerateGrowth_)
        {
            agglomerateCrystallization[momenti] +=
                agglomerateGrowthModel_->phaseSpaceConvection
                (
                    agglomerateOrders[momenti],
                    celli,
                    agglomerateQuadrature_
                );
        }
    }

    scalar crystallizationScale = 1.0;

    if (speciesCoupled_)
    {
        const label singleM3 = momentIndex(singleCrystalQuadrature_, 3);
        const label agglomerateM3 =
            momentIndex(agglomerateQuadrature_, 3);
        const scalar solidMassRate =
            rhop_.value()*shapeFactor_.value()
           *max
            (
                singleCrystallization[singleM3]
              + agglomerateCrystallization[agglomerateM3],
                scalar(0)
            );

        if (solidMassRate > SMALL)
        {
            const volScalarField& C =
                phi_.mesh().lookupObject<volScalarField>("C");
            crystallizationScale = min
            (
                scalar(1),
                0.9*max(C[celli], scalar(0))
               /(
                    max(phi_.mesh().time().deltaTValue(), scalar(SMALL))
                   *solidMassRate
                )
            );
        }
    }

    forAll(singleSources, momenti)
    {
        singleSources[momenti] =
            crystallizationScale*singleCrystallization[momenti];
    }
    forAll(agglomerateSources, momenti)
    {
        agglomerateSources[momenti] =
            crystallizationScale*agglomerateCrystallization[momenti];
    }

    updateCrossGrowthRate();

    const mappedPtrList<volScalarNode>& singleNodes =
        singleCrystalQuadrature_.nodes();
    const mappedPtrList<volScalarNode>& agglomerateNodes =
        agglomerateQuadrature_.nodes();

    forAll(singleAggregation, momenti)
    {
        const label order = singleOrders[momenti][0];
        scalar source = 0.0;

        if (ssEnabled_)
        {
            forAll(singleNodes, node1i)
            {
                const volScalarNode& node1 = singleNodes[node1i];
                const scalar L1 =
                    max(node1.abscissae()[0][celli], scalar(0));
                const scalar n1 = max
                (
                    node1.numberDensity
                    (
                        celli,
                        node1.weight()[celli],
                        L1
                    ),
                    scalar(0)
                );

                forAll(singleNodes, node2i)
                {
                    const volScalarNode& node2 = singleNodes[node2i];
                    const scalar L2 =
                        max(node2.abscissae()[0][celli], scalar(0));
                    const scalar n2 = max
                    (
                        node2.numberDensity
                        (
                            celli,
                            node2.weight()[celli],
                            L2
                        ),
                        scalar(0)
                    );

                    source -=
                        n1*n2
                       *aggregationKernelValue
                        (
                            *ssKernel_,
                            node1,
                            node2,
                            celli
                        )
                       *pow(L1, order);
                }
            }
        }

        if (saEnabled_)
        {
            forAll(singleNodes, node1i)
            {
                const volScalarNode& node1 = singleNodes[node1i];
                const scalar L1 =
                    max(node1.abscissae()[0][celli], scalar(0));
                const scalar n1 = max
                (
                    node1.numberDensity
                    (
                        celli,
                        node1.weight()[celli],
                        L1
                    ),
                    scalar(0)
                );

                forAll(agglomerateNodes, node2i)
                {
                    const volScalarNode& node2 =
                        agglomerateNodes[node2i];
                    const scalar L2 =
                        max(node2.abscissae()[0][celli], scalar(0));
                    const scalar n2 = max
                    (
                        node2.numberDensity
                        (
                            celli,
                            node2.weight()[celli],
                            L2
                        ),
                        scalar(0)
                    );

                    source -=
                        n1*n2
                       *aggregationKernelValue
                        (
                            *saKernel_,
                            node1,
                            node2,
                            celli
                        )
                       *pow(L1, order);
                }
            }
        }

        singleAggregation[momenti] = source;
        singleSources[momenti] += source;
    }

    forAll(agglomerateSources, momenti)
    {
        const label order = agglomerateOrders[momenti][0];
        scalar source = 0.0;

        if (order == 3)
        {
            source =
               -singleAggregation
                [
                    momentIndex(singleCrystalQuadrature_, 3)
                ];
        }
        else
        {
            if (ssEnabled_)
            {
                forAll(singleNodes, node1i)
                {
                    const volScalarNode& node1 = singleNodes[node1i];
                    const scalar L1 =
                        max(node1.abscissae()[0][celli], scalar(0));
                    const scalar n1 = max
                    (
                        node1.numberDensity
                        (
                            celli,
                            node1.weight()[celli],
                            L1
                        ),
                        scalar(0)
                    );

                    forAll(singleNodes, node2i)
                    {
                        const volScalarNode& node2 =
                            singleNodes[node2i];
                        const scalar L2 =
                            max(node2.abscissae()[0][celli], scalar(0));
                        const scalar n2 = max
                        (
                            node2.numberDensity
                            (
                                celli,
                                node2.weight()[celli],
                                L2
                            ),
                            scalar(0)
                        );
                        const scalar Lnew =
                            cbrt(pow3(L1) + pow3(L2));

                        source +=
                            0.5*n1*n2
                           *aggregationKernelValue
                            (
                                *ssKernel_,
                                node1,
                                node2,
                                celli
                            )
                           *pow(Lnew, order);
                    }
                }
            }

            if (saEnabled_)
            {
                forAll(singleNodes, node1i)
                {
                    const volScalarNode& node1 = singleNodes[node1i];
                    const scalar L1 =
                        max(node1.abscissae()[0][celli], scalar(0));
                    const scalar n1 = max
                    (
                        node1.numberDensity
                        (
                            celli,
                            node1.weight()[celli],
                            L1
                        ),
                        scalar(0)
                    );

                    forAll(agglomerateNodes, node2i)
                    {
                        const volScalarNode& node2 =
                            agglomerateNodes[node2i];
                        const scalar L2 =
                            max(node2.abscissae()[0][celli], scalar(0));
                        const scalar n2 = max
                        (
                            node2.numberDensity
                            (
                                celli,
                                node2.weight()[celli],
                                L2
                            ),
                            scalar(0)
                        );
                        const scalar Lnew =
                            cbrt(pow3(L1) + pow3(L2));

                        source +=
                            n1*n2
                           *aggregationKernelValue
                            (
                                *saKernel_,
                                node1,
                                node2,
                                celli
                            )
                           *(pow(Lnew, order) - pow(L2, order));
                    }
                }
            }

            if (aaEnabled_)
            {
                forAll(agglomerateNodes, node1i)
                {
                    const volScalarNode& node1 =
                        agglomerateNodes[node1i];
                    const scalar L1 =
                        max(node1.abscissae()[0][celli], scalar(0));
                    const scalar n1 = max
                    (
                        node1.numberDensity
                        (
                            celli,
                            node1.weight()[celli],
                            L1
                        ),
                        scalar(0)
                    );

                    forAll(agglomerateNodes, node2i)
                    {
                        const volScalarNode& node2 =
                            agglomerateNodes[node2i];
                        const scalar L2 =
                            max(node2.abscissae()[0][celli], scalar(0));
                        const scalar n2 = max
                        (
                            node2.numberDensity
                            (
                                celli,
                                node2.weight()[celli],
                                L2
                            ),
                            scalar(0)
                        );
                        const scalar Lnew =
                            cbrt(pow3(L1) + pow3(L2));

                        source +=
                            0.5*n1*n2
                           *aggregationKernelValue
                            (
                                *aaKernel_,
                                node1,
                                node2,
                                celli
                            )
                           *(
                                pow(Lnew, order)
                              - pow(L1, order)
                              - pow(L2, order)
                            );
                    }
                }
            }
        }

        agglomerateSources[momenti] += source;
    }

    if (agglomerateBreakup_)
    {
        addBreakupSources(celli, singleSources, agglomerateSources);
    }

    forAll(singleSources, momenti)
    {
        if (!std::isfinite(singleSources[momenti]))
        {
            FatalErrorInFunction
                << "Non-finite singleCrystal source for moment "
                << singleOrders[momenti] << " in cell " << celli
                << abort(FatalError);
        }
    }

    forAll(agglomerateSources, momenti)
    {
        if (!std::isfinite(agglomerateSources[momenti]))
        {
            FatalErrorInFunction
                << "Non-finite agglomerate source for moment "
                << agglomerateOrders[momenti] << " in cell " << celli
                << abort(FatalError);
        }
    }
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::addBreakupSources
(
    const label celli,
    scalarList& singleSources,
    scalarList& agglomerateSources
) const
{
    const auto& nodes = agglomerateQuadrature_.nodes();
    const auto& singleOrders = singleCrystalQuadrature_.momentOrders();
    const auto& agglomerateOrders = agglomerateQuadrature_.momentOrders();
    const scalar f = singleCrystalDaughterFraction_;
    scalar volumeTransfer = 0;

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        const scalar L = max(node.abscissae()[0][celli], scalar(0));
        const scalar number = node.numberDensity(celli, node.weight()[celli], L);
        if (number <= 0 || L <= 0)
        {
            continue;
        }
        const scalar events = number*agglomerateBreakupKernel_->Kb(L, celli);
        volumeTransfer += f*events*pow3(L);

        forAll(singleSources, momenti)
        {
            const label k = singleOrders[momenti][0];
            if (k != 3 && f > 0)
            {
                singleSources[momenti] += f*events
                    *agglomerateBreakupKernel_->daughterMoment(k, L);
            }
        }
        forAll(agglomerateSources, momenti)
        {
            const label k = agglomerateOrders[momenti][0];
            if (k != 3)
            {
                agglomerateSources[momenti] += events
                    *((1 - f)*agglomerateBreakupKernel_->daughterMoment(k, L)
                      - Foam::pow(L, k));
            }
        }
    }
    // A common transfer enforces total solid-volume conservation exactly.
    singleSources[momentIndex(singleCrystalQuadrature_, 3)] += volumeTransfer;
    agglomerateSources[momentIndex(agglomerateQuadrature_, 3)] -= volumeTransfer;
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::setCellMoments
(
    const label celli,
    const scalarList& singleValues,
    const scalarList& agglomerateValues
)
{
    volScalarMomentFieldSet& singleMoments =
        singleCrystalQuadrature_.moments();
    volScalarMomentFieldSet& agglomerateMoments =
        agglomerateQuadrature_.moments();

    forAll(singleValues, momenti)
    {
        singleMoments[momenti][celli] = singleValues[momenti];
    }
    forAll(agglomerateValues, momenti)
    {
        agglomerateMoments[momenti][celli] =
            agglomerateValues[momenti];
    }
}


bool Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::updateCellQuadratures
(
    const label celli,
    const bool fatal
)
{
    // The generic non-fatal inversion rejects m0=0. An exactly empty
    // population is valid here and must not block the other population.
    const auto invert = [celli, fatal](scalarQuadratureApproximation& q)
    {
        bool empty = true;
        forAll(q.moments(), momenti)
        {
            empty = empty && q.moments()[momenti][celli] == 0;
        }
        if (empty)
        {
            // updateLocalQuadrature always invokes the generic non-fatal
            // inverter, so reset the empty node representation explicitly.
            forAll(q.nodes(), nodei)
            {
                q.nodes()[nodei].weight()[celli] = 0;
                q.nodes()[nodei].abscissae()[0][celli] = 0;
            }
            return true;
        }
        return q.updateLocalQuadrature(celli, fatal);
    };
    const bool singleRealizable = invert(singleCrystalQuadrature_);
    const bool agglomerateRealizable = invert(agglomerateQuadrature_);
    const bool bothRealizable =
        singleRealizable && agglomerateRealizable;

    if (bothRealizable)
    {
        singleCrystalQuadrature_.updateLocalMoments(celli);
        agglomerateQuadrature_.updateLocalMoments(celli);
    }

    return bothRealizable;
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::solveCoupledSources()
{
    if
    (
        !solveSources_
     ||
        (
            !nucleation_
         && !singleGrowth_
         && !agglomerateGrowth_
         && !ssEnabled_
         && !saEnabled_
         && !aaEnabled_
         && !agglomerateBreakup_
        )
    )
    {
        return;
    }

    volScalarMomentFieldSet& singleMoments =
        singleCrystalQuadrature_.moments();
    volScalarMomentFieldSet& agglomerateMoments =
        agglomerateQuadrature_.moments();

    const label nSingle = singleMoments.size();
    const label nAgglomerate = agglomerateMoments.size();
    const scalar globalDt = phi_.mesh().time().deltaTValue();

    Info<< "Solving two-population source terms with a common realizable "
        << "integrator." << endl;

    forAll(singleMoments[0], celli)
    {
        singleCrystalQuadrature_.updateLocalQuadrature(celli);
        agglomerateQuadrature_.updateLocalQuadrature(celli);

        scalarList oldSingle(nSingle, 0.0);
        scalarList oldAgglomerate(nAgglomerate, 0.0);

        forAll(oldSingle, momenti)
        {
            oldSingle[momenti] = singleMoments[momenti][celli];
        }
        forAll(oldAgglomerate, momenti)
        {
            oldAgglomerate[momenti] =
                agglomerateMoments[momenti][celli];
        }

        if (!solveOde_)
        {
            scalarList singleSources;
            scalarList agglomerateSources;
            calculateSources
            (
                celli,
                singleSources,
                agglomerateSources
            );

            scalarList newSingle(oldSingle);
            scalarList newAgglomerate(oldAgglomerate);

            forAll(newSingle, momenti)
            {
                newSingle[momenti] +=
                    globalDt*singleSources[momenti];
            }
            forAll(newAgglomerate, momenti)
            {
                newAgglomerate[momenti] +=
                    globalDt*agglomerateSources[momenti];
            }

            setCellMoments(celli, newSingle, newAgglomerate);
            updateCellQuadratures(celli, true);
            continue;
        }

        scalar localT = 0.0;
        scalar localDt = min(localDt_[celli], globalDt);

        while
        (
            localT
          < globalDt - SMALL*max(globalDt, scalar(1))
        )
        {
            localDt = min(localDt, globalDt - localT);

            bool accepted = false;

            while (!accepted)
            {
                scalarList k1Single(nSingle, 0.0);
                scalarList k2Single(nSingle, 0.0);
                scalarList k3Single(nSingle, 0.0);
                scalarList k1Agglomerate(nAgglomerate, 0.0);
                scalarList k2Agglomerate(nAgglomerate, 0.0);
                scalarList k3Agglomerate(nAgglomerate, 0.0);
                scalarList sourceSingle;
                scalarList sourceAgglomerate;
                scalarList trialSingle(oldSingle);
                scalarList trialAgglomerate(oldAgglomerate);

                calculateSources
                (
                    celli,
                    sourceSingle,
                    sourceAgglomerate
                );

                forAll(k1Single, momenti)
                {
                    k1Single[momenti] =
                        localDt*sourceSingle[momenti];
                    trialSingle[momenti] =
                        oldSingle[momenti] + k1Single[momenti];
                }
                forAll(k1Agglomerate, momenti)
                {
                    k1Agglomerate[momenti] =
                        localDt*sourceAgglomerate[momenti];
                    trialAgglomerate[momenti] =
                        oldAgglomerate[momenti]
                      + k1Agglomerate[momenti];
                }

                setCellMoments
                (
                    celli,
                    trialSingle,
                    trialAgglomerate
                );

                bool realizable =
                    updateCellQuadratures(celli, false);

                if (realizable)
                {
                    calculateSources
                    (
                        celli,
                        sourceSingle,
                        sourceAgglomerate
                    );

                    forAll(k2Single, momenti)
                    {
                        k2Single[momenti] =
                            localDt*sourceSingle[momenti];
                        trialSingle[momenti] =
                            oldSingle[momenti]
                          + 0.25
                           *(
                                k1Single[momenti]
                              + k2Single[momenti]
                            );
                    }
                    forAll(k2Agglomerate, momenti)
                    {
                        k2Agglomerate[momenti] =
                            localDt*sourceAgglomerate[momenti];
                        trialAgglomerate[momenti] =
                            oldAgglomerate[momenti]
                          + 0.25
                           *(
                                k1Agglomerate[momenti]
                              + k2Agglomerate[momenti]
                            );
                    }

                    setCellMoments
                    (
                        celli,
                        trialSingle,
                        trialAgglomerate
                    );
                    realizable = updateCellQuadratures(celli, false);
                }

                if (realizable)
                {
                    calculateSources
                    (
                        celli,
                        sourceSingle,
                        sourceAgglomerate
                    );

                    forAll(k3Single, momenti)
                    {
                        k3Single[momenti] =
                            localDt*sourceSingle[momenti];
                        trialSingle[momenti] =
                            oldSingle[momenti]
                          + (
                                k1Single[momenti]
                              + k2Single[momenti]
                              + 4.0*k3Single[momenti]
                            )/6.0;
                    }
                    forAll(k3Agglomerate, momenti)
                    {
                        k3Agglomerate[momenti] =
                            localDt*sourceAgglomerate[momenti];
                        trialAgglomerate[momenti] =
                            oldAgglomerate[momenti]
                          + (
                                k1Agglomerate[momenti]
                              + k2Agglomerate[momenti]
                              + 4.0*k3Agglomerate[momenti]
                            )/6.0;
                    }

                    setCellMoments
                    (
                        celli,
                        trialSingle,
                        trialAgglomerate
                    );
                    realizable = updateCellQuadratures(celli, false);
                }

                if (!realizable)
                {
                    if (localDtAdjustments_ == 0)
                    {
                        Info<< "A coupled moment stage was not realizable; "
                            << "reducing the common local time step." << endl;
                    }
                    ++localDtAdjustments_;
                    setCellMoments
                    (
                        celli,
                        oldSingle,
                        oldAgglomerate
                    );
                    updateCellQuadratures(celli, true);
                    localDt *= 0.5;

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minLocalDt while preserving both "
                            << "population moment sets in cell " << celli
                            << abort(FatalError);
                    }
                    continue;
                }

                scalar errorSum = 0.0;
                label nError = 0;

                forAll(k3Single, momenti)
                {
                    const scalar diff =
                        (
                            2.0*k3Single[momenti]
                          - k1Single[momenti]
                          - k2Single[momenti]
                        )/3.0;
                    const scalar scale =
                        ATol_
                      + RTol_*max
                        (
                            mag(singleMoments[momenti][celli]),
                            mag(oldSingle[momenti])
                        );
                    errorSum += sqr(diff/scale);
                    ++nError;
                }

                forAll(k3Agglomerate, momenti)
                {
                    const scalar diff =
                        (
                            2.0*k3Agglomerate[momenti]
                          - k1Agglomerate[momenti]
                          - k2Agglomerate[momenti]
                        )/3.0;
                    const scalar scale =
                        ATol_
                      + RTol_*max
                        (
                            mag(agglomerateMoments[momenti][celli]),
                            mag(oldAgglomerate[momenti])
                        );
                    errorSum += sqr(diff/scale);
                    ++nError;
                }

                const scalar error =
                    sqrt(errorSum/max(nError, label(1)));

                if (error <= 1.0)
                {
                    accepted = true;
                    localT += localDt;

                    forAll(oldSingle, momenti)
                    {
                        oldSingle[momenti] =
                            singleMoments[momenti][celli];
                    }
                    forAll(oldAgglomerate, momenti)
                    {
                        oldAgglomerate[momenti] =
                            agglomerateMoments[momenti][celli];
                    }

                    const scalar factor =
                        error < SMALL
                      ? facMax_
                      : min
                        (
                            facMax_,
                            max
                            (
                                facMin_,
                                fac_/pow(error, 1.0/3.0)
                            )
                        );
                    localDt *= factor;
                    localDt_[celli] = localDt;
                }
                else
                {
                    setCellMoments
                    (
                        celli,
                        oldSingle,
                        oldAgglomerate
                    );
                    updateCellQuadratures(celli, true);
                    localDt *= min
                    (
                        scalar(1),
                        max
                        (
                            facMin_,
                            fac_/pow(error, 1.0/3.0)
                        )
                    );

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minLocalDt while satisfying the "
                            << "coupled source error tolerance in cell "
                            << celli << abort(FatalError);
                    }
                }
            }
        }
    }

    forAll(singleMoments, momenti)
    {
        singleMoments[momenti].correctBoundaryConditions();
    }
    forAll(agglomerateMoments, momenti)
    {
        agglomerateMoments[momenti].correctBoundaryConditions();
    }

    singleCrystalQuadrature_.updateBoundaryQuadrature();
    agglomerateQuadrature_.updateBoundaryQuadrature();
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::updateStatistics()
{
    updateTotalM2();

    const volScalarMomentFieldSet& singleMoments =
        singleCrystalQuadrature_.moments();
    const volScalarMomentFieldSet& agglomerateMoments =
        agglomerateQuadrature_.moments();

    const volScalarMoment& s0 = singleMoments(labelList(1, 0));
    const volScalarMoment& s1 = singleMoments(labelList(1, 1));
    const volScalarMoment& s2 = singleMoments(labelList(1, 2));
    const volScalarMoment& s3 = singleMoments(labelList(1, 3));
    const volScalarMoment& a0 = agglomerateMoments(labelList(1, 0));
    const volScalarMoment& a1 = agglomerateMoments(labelList(1, 1));
    const volScalarMoment& a2 = agglomerateMoments(labelList(1, 2));
    const volScalarMoment& a3 = agglomerateMoments(labelList(1, 3));

    forAll(singleL10_, celli)
    {
        singleL10_[celli] =
            s1[celli]/max(s0[celli], scalar(SMALL));
        singleL32_[celli] =
            s3[celli]/max(s2[celli], scalar(SMALL));
        agglomerateL10_[celli] =
            a1[celli]/max(a0[celli], scalar(SMALL));
        agglomerateL32_[celli] =
            a3[celli]/max(a2[celli], scalar(SMALL));
        totalL10_[celli] =
            (s1[celli] + a1[celli])
           /max(s0[celli] + a0[celli], scalar(SMALL));
        totalSolidMass_[celli] =
            rhop_.value()*shapeFactor_.value()
           *max(s3[celli] + a3[celli], scalar(0));
    }

    singleL10_.correctBoundaryConditions();
    singleL32_.correctBoundaryConditions();
    agglomerateL10_.correctBoundaryConditions();
    agglomerateL32_.correctBoundaryConditions();
    totalL10_.correctBoundaryConditions();
    totalSolidMass_.correctBoundaryConditions();
    updateCrossGrowthRate();
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::updateObservedMoments()
{
    if (!observedMoment_)
    {
        return;
    }

    forAll(observedMoments_, momenti)
    {
        observedMoments_[momenti] = dimensionedScalar
        (
            "zero",
            observedMoments_[momenti].dimensions(),
            0.0
        );
    }

    const scalar logL50 = Foam::log(detectionL50_.value());

    const scalarQuadratureApproximation* quadratures[2] =
    {
        &singleCrystalQuadrature_,
        &agglomerateQuadrature_
    };

    for (label populationi = 0; populationi < 2; ++populationi)
    {
        const mappedPtrList<volScalarNode>& nodes =
            quadratures[populationi]->nodes();

        forAll(nodes, nodei)
        {
            const volScalarField& weight = nodes[nodei].weight();
            const volScalarField& length = nodes[nodei].abscissae()[0];

            forAll(weight, celli)
            {
                const scalar Li = max(length[celli], VSMALL);
                const scalar logDistance =
                    (Foam::log(Li) - logL50)/detectionLogWidth_;
                const scalar boundedDistance =
                    max(min(logDistance, scalar(60)), scalar(-60));
                const scalar detectionProbability =
                    1.0/(1.0 + Foam::exp(-boundedDistance));

                forAll(observedMoments_, momenti)
                {
                    observedMoments_[momenti][celli] +=
                        weight[celli]
                       *detectionProbability
                       *Foam::pow
                        (
                            Li,
                            observedMomentOrders_[momenti]
                        );
                }
            }
        }
    }

    label m2Index = -1;
    label m3Index = -1;

    forAll(observedMomentOrders_, momenti)
    {
        if (observedMomentOrders_[momenti] == 2)
        {
            m2Index = momenti;
        }
        else if (observedMomentOrders_[momenti] == 3)
        {
            m3Index = momenti;
        }
    }

    forAll(observedD32_(), celli)
    {
        observedD32_()[celli] =
            m2Index >= 0 && m3Index >= 0
          ? observedMoments_[m3Index][celli]
           /(observedMoments_[m2Index][celli] + VSMALL)
          : 0.0;
    }

    forAll(observedMoments_, momenti)
    {
        observedMoments_[momenti].correctBoundaryConditions();
    }
    observedD32_().correctBoundaryConditions();
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::calculateSpeciesTransfer
(
    const scalarField& totalM3Before
)
{
    SYact_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);

    if (!speciesCoupled_ || !solveSources_)
    {
        SYact_.correctBoundaryConditions();
        return;
    }

    const volScalarMoment& singleM3 =
        singleCrystalQuadrature_.moments()(labelList(1, 3));
    const volScalarMoment& agglomerateM3 =
        agglomerateQuadrature_.moments()(labelList(1, 3));
    const scalar invDt =
        1.0/max(phi_.mesh().time().deltaTValue(), scalar(SMALL));
    const scalar factor = rhop_.value()*shapeFactor_.value();

    forAll(SYact_, celli)
    {
        SYact_[celli] =
            factor
           *(
                singleM3[celli]
              + agglomerateM3[celli]
              - totalM3Before[celli]
            )*invDt;
    }

    SYact_.correctBoundaryConditions();
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::realizableCo() const
{
    return min
    (
        singleCrystalAdvection_->realizableCo(),
        agglomerateAdvection_->realizableCo()
    );
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::CoNum() const
{
    return 0.0;
}


void Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::solve()
{
    preUpdate();

    singleCrystalAdvection_->update();
    agglomerateAdvection_->update();

    fv::EulerDdtScheme<scalar> momentDdt(phi_.mesh());

    volScalarMomentFieldSet& singleMoments =
        singleCrystalQuadrature_.moments();
    volScalarMomentFieldSet& agglomerateMoments =
        agglomerateQuadrature_.moments();

    forAll(singleMoments, momenti)
    {
        volScalarMoment& moment = singleMoments[momenti];
        fvScalarMatrix equation
        (
            momentDdt.fvmDdt(moment)
          + singleCrystalAdvection_->divMoments()[momenti]
         == singleDiffusionModel_->momentDiff(moment)
        );

        if (singleDiffusionModel_->type() != "none")
        {
            equation.relax();
        }

        equation.solve();
        moment.correctBoundaryConditions();
    }

    forAll(agglomerateMoments, momenti)
    {
        volScalarMoment& moment = agglomerateMoments[momenti];
        fvScalarMatrix equation
        (
            momentDdt.fvmDdt(moment)
          + agglomerateAdvection_->divMoments()[momenti]
         == agglomerateDiffusionModel_->momentDiff(moment)
        );

        if (agglomerateDiffusionModel_->type() != "none")
        {
            equation.relax();
        }

        equation.solve();
        moment.correctBoundaryConditions();
    }

    singleCrystalQuadrature_.updateQuadrature();
    agglomerateQuadrature_.updateQuadrature();

    const scalarField totalM3Before
    (
        singleMoments(labelList(1, 3)).primitiveField()
      + agglomerateMoments(labelList(1, 3)).primitiveField()
    );

    solveCoupledSources();
    updateStatistics();
    updateObservedMoments();
    calculateSpeciesTransfer(totalM3Before);
}


bool Foam::PDFTransportModels::populationBalanceModels::
twoPopulationBalanceSystem::readIfModified()
{
    const dictionary& dict =
        populationBalanceProperties_.subDict(type() + "Coeffs");
    const dictionary& odeDict = dict.subDict("odeCoeffs");

    solveSources_ = odeDict.lookupOrDefault<Switch>("solveSources", true);
    solveOde_ = odeDict.lookupOrDefault<Switch>("solveOde", true);
    odeDict.lookup("ATol") >> ATol_;
    odeDict.lookup("RTol") >> RTol_;
    odeDict.lookup("fac") >> fac_;
    odeDict.lookup("facMin") >> facMin_;
    odeDict.lookup("facMax") >> facMax_;
    odeDict.lookup("minLocalDt") >> minLocalDt_;

    return true;
}

// ************************************************************************* //

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

#include "crysDualImpurityPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"
#include "EulerDdtScheme.H"
#include "zeroGradientFvPatchField.H"
#include "mathematicalConstants.H"
#include <cmath>

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(crysDualImpurityPopulationBalance, 0);

    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        crysDualImpurityPopulationBalance,
        dictionary
    );
}
}
}


Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::crysDualImpurityPopulationBalance
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
        List<supportType>(3, supportType::RPlus)
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
    nucleation_(dict.lookupOrDefault("nucleation", true)),
    growth_(dict.lookupOrDefault("growth", true)),
    speciesCoupled_(dict.lookupOrDefault("speciesCoupled", true)),
    nucleationModel_(),
    diffusionModel_
    (
        populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    T_
    (
        phi.mesh().lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("TField", "T")
        )
    ),
    C_
    (
        phi.mesh().lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("soluteField", "C")
        )
    ),
    Ci1_
    (
        phi.mesh().lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("impurity1Field", "Ci1")
        )
    ),
    Ci2_
    (
        phi.mesh().lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("impurity2Field", "Ci2")
        )
    ),
    sigma_
    (
        phi.mesh().lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("supersaturationField", "sigma")
        )
    ),
    rhop_
    (
        dimensionedScalar::getOrDefault
        (
            "rhop",
            dict,
            dimDensity,
            2338.0
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
            1.0e-4
        )
    ),
    mRef_
    (
        "mRef",
        dimMass,
        rhop_.value()*shapeFactor_.value()*pow3(sizeRef_.value())
    ),
    rhoSolvent_
    (
        dimensionedScalar::getOrDefault
        (
            "rhoSolvent",
            dict,
            dimDensity,
            1000.0
        )
    ),
    kg1_("kg1", dimLength/dimTime, dict.subDict("growthCoeffs")),
    kg2_("kg2", dimLength/dimTime, dict.subDict("growthCoeffs")),
    g1_(readScalar(dict.subDict("growthCoeffs").lookup("g1"))),
    g2_(readScalar(dict.subDict("growthCoeffs").lookup("g2"))),
    growthScale_
    (
        dict.subDict("reductionCoeffs").lookupOrDefault<scalar>
        (
            "growthScale",
            1.0
        )
    ),
    kads01_
    (
        readScalar(dict.subDict("adsorptionCoeffs").lookup("kads01"))
    ),
    kdes01_
    (
        readScalar(dict.subDict("adsorptionCoeffs").lookup("kdes01"))
    ),
    kads02_
    (
        readScalar(dict.subDict("adsorptionCoeffs").lookup("kads02"))
    ),
    kdes02_
    (
        readScalar(dict.subDict("adsorptionCoeffs").lookup("kdes02"))
    ),
    deltaGads1_
    (
        "deltaGads1",
        dimEnergy/dimMoles,
        dict.subDict("adsorptionCoeffs")
    ),
    deltaGdes1_
    (
        "deltaGdes1",
        dimEnergy/dimMoles,
        dict.subDict("adsorptionCoeffs")
    ),
    deltaGads2_
    (
        "deltaGads2",
        dimEnergy/dimMoles,
        dict.subDict("adsorptionCoeffs")
    ),
    deltaGdes2_
    (
        "deltaGdes2",
        dimEnergy/dimMoles,
        dict.subDict("adsorptionCoeffs")
    ),
    gasConstant_
    (
        dimensionedScalar::getOrDefault
        (
            "R",
            dict.subDict("adsorptionCoeffs"),
            dimEnergy/dimMoles/dimTemperature,
            8.31446261815324
        )
    ),
    beta1_
    (
        "beta1",
        dimTemperature,
        dict.subDict("adsorptionCoeffs")
    ),
    beta2_
    (
        "beta2",
        dimTemperature,
        dict.subDict("adsorptionCoeffs")
    ),
    faceEffectScale1_
    (
        dict.subDict("reductionCoeffs").lookupOrDefault<scalar>
        (
            "faceEffectScale1",
            1.0
        )
    ),
    faceEffectScale2_
    (
        dict.subDict("reductionCoeffs").lookupOrDefault<scalar>
        (
            "faceEffectScale2",
            1.0
        )
    ),
    molarMassKDP_
    (
        "molarMassKDP",
        dimMass/dimMoles,
        dict.subDict("partitionCoeffs")
    ),
    molarMassCGM1_
    (
        "molarMassCGM1",
        dimMass/dimMoles,
        dict.subDict("partitionCoeffs")
    ),
    molarMassCGM2_
    (
        "molarMassCGM2",
        dimMass/dimMoles,
        dict.subDict("partitionCoeffs")
    ),
    Gmin1_("Gmin1", dimLength/dimTime, dict.subDict("partitionCoeffs")),
    Gmin2_("Gmin2", dimLength/dimTime, dict.subDict("partitionCoeffs")),
    km01_("km01", dimLength/dimTime, dict.subDict("partitionCoeffs")),
    km02_("km02", dimLength/dimTime, dict.subDict("partitionCoeffs")),
    Ke1_(readScalar(dict.subDict("partitionCoeffs").lookup("Ke1"))),
    Ke2_(readScalar(dict.subDict("partitionCoeffs").lookup("Ke2"))),
    partitionScale1_
    (
        dict.subDict("reductionCoeffs").lookupOrDefault<scalar>
        (
            "partitionScale1",
            1.0
        )
    ),
    partitionScale2_
    (
        dict.subDict("reductionCoeffs").lookupOrDefault<scalar>
        (
            "partitionScale2",
            1.0
        )
    ),
    z1Nucleation_
    (
        dimensionedScalar::getOrDefault
        (
            "z1Nucleation",
            dict,
            dimless,
            0.0
        )
    ),
    z2Nucleation_
    (
        dimensionedScalar::getOrDefault
        (
            "z2Nucleation",
            dict,
            dimless,
            0.0
        )
    ),
    L10_
    (
        IOobject
        (
            "L10.crysDualImpurity",
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
            "L32.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    meanZ1_
    (
        IOobject
        (
            "meanZ1.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    meanZ2_
    (
        IOobject
        (
            "meanZ2.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    solidCGM1Ppm_
    (
        IOobject
        (
            "solidCGM1Ppm.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    solidCGM2Ppm_
    (
        IOobject
        (
            "solidCGM2Ppm.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    totalImpurityPpm_
    (
        IOobject
        (
            "totalImpurityPpm.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    crystalPurityPercent_
    (
        IOobject
        (
            "crystalPurityPercent.crysDualImpurity",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimless, 100.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    faceFactor1_
    (
        IOobject
        (
            "borsosFaceFactor1",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("one", dimless, 1.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    faceFactor2_
    (
        IOobject
        (
            "borsosFaceFactor2",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("one", dimless, 1.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    faceGrowth1_
    (
        IOobject
        (
            "borsosFaceGrowth1",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    faceGrowth2_
    (
        IOobject
        (
            "borsosFaceGrowth2",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    equivalentGrowth_
    (
        IOobject
        (
            "equivalentGrowthRate",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0),
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
    SI1act_
    (
        IOobject
        (
            "SI1act",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    SI2act_
    (
        IOobject
        (
            "SI2act",
            phi.mesh().time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    crystallizationSourceScale_(phi.mesh().nCells(), 1.0)
{
    if
    (
        dict.lookupOrDefault("aggregation", false)
     || dict.lookupOrDefault("breakup", false)
    )
    {
        FatalIOErrorInFunction(dict)
            << "The reduced (L,z1,z2) model supports nucleation and growth "
            << "only; aggregation and breakup are not implemented."
            << exit(FatalIOError);
    }

    if (nucleation_)
    {
        nucleationModel_ = populationBalanceSubModels::nucleationModel::New
        (
            dict.subDict("nucleationModel"),
            phi.mesh()
        );
    }

    validateConfiguration();
    updateKineticDiagnostics();
    updateStatistics();

    const scalarField initialM300
    (
        quadrature_.moments()(labelList({3, 0, 0})).primitiveField()
    );
    const scalarField initialM010
    (
        quadrature_.moments()(labelList({0, 1, 0})).primitiveField()
    );
    const scalarField initialM001
    (
        quadrature_.moments()(labelList({0, 0, 1})).primitiveField()
    );
    calcSpeciesTransfer(initialM300, initialM010, initialM001);

    Info<< "Reduced dual-impurity crystal coordinates: (L,z1,z2), mRef="
        << mRef_.value() << " kg; volume-equivalent growth projection is "
        << "(6/pi)^(1/3)*(G1 + 2*G2)/3." << endl;
}


Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::~crysDualImpurityPopulationBalance()
{}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::validateConfiguration() const
{
    if (quadrature_.nDimensions() != 3)
    {
        FatalErrorInFunction
            << "crysDualImpurityPopulationBalance requires exactly three "
            << "internal coordinates (L,z1,z2). Found "
            << quadrature_.nDimensions() << abort(FatalError);
    }

    const labelList& nNodes = quadrature_.nNodes();
    if
    (
        nNodes.size() != 3
     || nNodes[0] < 2
     || nNodes[1] != 1
     || nNodes[2] != 1
    )
    {
        FatalErrorInFunction
            << "The endpoint model requires nNodes=[nL,1,1], nL>=2. Found "
            << nNodes << abort(FatalError);
    }

    const label nL = nNodes[0];
    labelListList expected(4*nL);
    label index = 0;
    for (label i = 0; i < 2*nL; ++i)
    {
        expected[index++] = labelList({i, 0, 0});
    }
    for (label i = 0; i < nL; ++i)
    {
        expected[index++] = labelList({i, 1, 0});
    }
    for (label i = 0; i < nL; ++i)
    {
        expected[index++] = labelList({i, 0, 1});
    }

    const labelListList& orders = quadrature_.momentOrders();
    if (orders.size() != expected.size())
    {
        FatalErrorInFunction
            << "nNodes=" << nNodes << " requires moments " << expected
            << ", but found " << orders << abort(FatalError);
    }

    forAll(expected, requiredi)
    {
        bool found = false;
        forAll(orders, orderi)
        {
            found = found || orders[orderi] == expected[requiredi];
        }
        if (!found)
        {
            FatalErrorInFunction
                << "Missing required moment " << expected[requiredi]
                << " in " << orders << abort(FatalError);
        }
    }

    const labelListList& nodeIndexes = quadrature_.nodeIndexes();
    if (nodeIndexes.size() != nL)
    {
        FatalErrorInFunction
            << "nNodes=" << nNodes << " requires " << nL
            << " nodes, found " << nodeIndexes << abort(FatalError);
    }
    for (label i = 0; i < nL; ++i)
    {
        const labelList required({i, 0, 0});
        bool found = false;
        forAll(nodeIndexes, nodei)
        {
            found = found || nodeIndexes[nodei] == required;
        }
        if (!found)
        {
            FatalErrorInFunction
                << "Missing required node " << required << abort(FatalError);
        }
    }

    const dictionary& fieldInversion =
        quadrature_.subDict("basicScalarMomentInversion");
    const dictionary& advection = quadrature_.subDict("momentAdvection");
    const dictionary& faceInversion =
        advection.subDict("basicScalarMomentInversion");
    const wordList fieldSupports(fieldInversion.lookup("supports"));
    const wordList faceSupports(faceInversion.lookup("supports"));

    if
    (
        word(quadrature_.lookup("fieldMomentInversion"))
            != "basicScalarFieldMomentInversion"
     || word(fieldInversion.lookup("type")) != "conditional"
     || word(advection.lookup("univariateMomentAdvection"))
            != "multivariateFirstOrderKineticScalar"
     || word(faceInversion.lookup("type")) != "conditional"
     || fieldSupports != wordList(3, "RPlus")
     || faceSupports != wordList(3, "RPlus")
    )
    {
        FatalErrorInFunction
            << "The dual-impurity model requires conditional CQMOM with "
            << "supports (RPlus RPlus RPlus) for cell and face inversion."
            << abort(FatalError);
    }

    const volScalarMomentFieldSet& moments = quadrature_.moments();
    const volScalarMoment& m000 = moments(labelList({0, 0, 0}));
    if
    (
        m000.dimensions() != dimless/dimVolume
     || moments(labelList({1, 0, 0})).dimensions()/m000.dimensions()
            != dimLength
     || moments(labelList({0, 1, 0})).dimensions()/m000.dimensions()
            != dimless
     || moments(labelList({0, 0, 1})).dimensions()/m000.dimensions()
            != dimless
    )
    {
        FatalErrorInFunction
            << "Internal-coordinate dimensions must be (length, "
            << "dimensionless, dimensionless)." << abort(FatalError);
    }

    if
    (
        rhop_.value() <= 0
     || shapeFactor_.value() <= 0
     || sizeRef_.value() <= 0
     || rhoSolvent_.value() <= 0
     || kg1_.value() < 0
     || kg2_.value() < 0
     || g1_ < 0
     || g2_ < 0
     || growthScale_ < 0
     || kads01_ < 0
     || kads02_ < 0
     || kdes01_ <= 0
     || kdes02_ <= 0
     || gasConstant_.value() <= 0
     || beta1_.value() < 0
     || beta2_.value() < 0
     || faceEffectScale1_ < 0
     || faceEffectScale2_ < 0
     || molarMassKDP_.value() <= 0
     || molarMassCGM1_.value() <= 0
     || molarMassCGM2_.value() <= 0
     || Gmin1_.value() < 0
     || Gmin2_.value() < 0
     || km01_.value() <= 0
     || km02_.value() <= 0
     || Ke1_ < 0 || Ke1_ > 1
     || Ke2_ < 0 || Ke2_ > 1
     || partitionScale1_ < 0
     || partitionScale2_ < 0
     || z1Nucleation_.value() < 0
     || z2Nucleation_.value() < 0
    )
    {
        FatalErrorInFunction
            << "Invalid dual-impurity physical coefficient. Positive rate, "
            << "density, mass and scale values are required; Ke must be in "
            << "[0,1]." << abort(FatalError);
    }
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::massTransferCoefficient
(
    const scalar growthRate,
    const scalar km0
) const
{
    if (growthRate <= ROOTSMALL*km0)
    {
        return km0;
    }
    return growthRate/(-std::expm1(-growthRate/km0));
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::kineticState
(
    const label celli,
    scalar& p1,
    scalar& p2,
    scalar& G1,
    scalar& G2,
    scalar& GL,
    scalar& q1,
    scalar& q2
) const
{
    p1 = 1.0;
    p2 = 1.0;
    G1 = 0.0;
    G2 = 0.0;
    GL = 0.0;
    q1 = 0.0;
    q2 = 0.0;

    const scalar temperature = T_[celli];
    const scalar sigma = max(sigma_[celli], scalar(0));

    if (temperature <= 0 || sigma <= 1.0e-12)
    {
        return;
    }

    const scalar cgm1 = max(Ci1_[celli], scalar(0))/rhoSolvent_.value();
    const scalar cgm2 = max(Ci2_[celli], scalar(0))/rhoSolvent_.value();

    const scalar K1 =
        kads01_/kdes01_
       *std::exp
        (
            (deltaGdes1_.value() - deltaGads1_.value())
           /(gasConstant_.value()*temperature)
        );
    const scalar K2 =
        kads02_/kdes02_
       *std::exp
        (
            (deltaGdes2_.value() - deltaGads2_.value())
           /(gasConstant_.value()*temperature)
        );

    const scalar theta1 = K1*cgm1/(1.0 + K1*cgm1);
    const scalar theta2 = K2*cgm2/(1.0 + K2*cgm2);
    const scalar alpha1 =
        faceEffectScale1_*beta1_.value()/(temperature*sigma);
    const scalar alpha2 =
        faceEffectScale2_*beta2_.value()/(temperature*sigma);

    // A negative paper impurity factor represents complete growth arrest in
    // the reduced positive-growth model. Clamp only at the physical bounds.
    p1 = min(max(1.0 - alpha1*theta1, scalar(0)), scalar(1));
    p2 = min(max(1.0 - alpha2*theta2, scalar(0)), scalar(1));

    G1 = growthScale_*kg1_.value()*std::pow(sigma, g1_)*p1;
    G2 = growthScale_*kg2_.value()*std::pow(sigma, g2_)*p2;
    GL =
        std::pow(6.0/constant::mathematical::pi, 1.0/3.0)
       *(G1 + 2.0*G2)/3.0;

    if
    (
        !std::isfinite(K1) || !std::isfinite(K2)
     || !std::isfinite(G1) || !std::isfinite(G2) || !std::isfinite(GL)
    )
    {
        FatalErrorInFunction
            << "Non-finite Borsos kinetic state in cell " << celli
            << ": K1=" << K1 << ", K2=" << K2 << ", G1=" << G1
            << ", G2=" << G2 << ", GL=" << GL << exit(FatalError);
    }

    if (GL <= SMALL)
    {
        return;
    }

    const scalar km1 = massTransferCoefficient(G1, km01_.value());
    const scalar km2 = massTransferCoefficient(G2, km02_.value());
    const scalar kmin1 =
        massTransferCoefficient(Gmin1_.value(), km01_.value());
    const scalar kmin2 =
        massTransferCoefficient(Gmin2_.value(), km02_.value());

    const scalar exponent1 = G1 > SMALL
        ? Gmin1_.value()*km1/(G1*kmin1)
        : GREAT;
    const scalar exponent2 = G2 > SMALL
        ? Gmin2_.value()*km2/(G2*kmin2)
        : GREAT;
    const scalar base1 = max(1.0 - Ke1_, scalar(0));
    const scalar base2 = max(1.0 - Ke2_, scalar(0));
    const scalar Kd1 =
        (base1 <= SMALL || exponent1 >= 1000.0)
      ? 1.0
      : 1.0 - std::pow(base1, exponent1);
    const scalar Kd2 =
        (base2 <= SMALL || exponent2 >= 1000.0)
      ? 1.0
      : 1.0 - std::pow(base2, exponent2);

    const scalar soluteMoles =
        max(C_[celli], scalar(0))/molarMassKDP_.value();
    const scalar cgm1Moles =
        max(Ci1_[celli], scalar(0))/molarMassCGM1_.value();
    const scalar cgm2Moles =
        max(Ci2_[celli], scalar(0))/molarMassCGM2_.value();
    const scalar moleTotal = soluteMoles + cgm1Moles + cgm2Moles;

    if (moleTotal <= SMALL)
    {
        return;
    }

    scalar chi1 = partitionScale1_*Kd1*cgm1Moles/moleTotal;
    scalar chi2 = partitionScale2_*Kd2*cgm2Moles/moleTotal;
    const scalar chiSum = chi1 + chi2;

    if (chiSum >= 1.0 - ROOTSMALL)
    {
        const scalar scale = (1.0 - ROOTSMALL)/max(chiSum, scalar(SMALL));
        chi1 *= scale;
        chi2 *= scale;
    }

    const scalar denominator = max(1.0 - chi1 - chi2, scalar(ROOTSMALL));
    q1 =
        chi1/denominator
       *molarMassCGM1_.value()/molarMassKDP_.value();
    q2 =
        chi2/denominator
       *molarMassCGM2_.value()/molarMassKDP_.value();
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::updateKineticDiagnostics()
{
    forAll(faceFactor1_, celli)
    {
        scalar p1, p2, G1, G2, GL, q1, q2;
        kineticState(celli, p1, p2, G1, G2, GL, q1, q2);
        faceFactor1_[celli] = p1;
        faceFactor2_[celli] = p2;
        faceGrowth1_[celli] = G1;
        faceGrowth2_[celli] = G2;
        equivalentGrowth_[celli] = GL;
    }
    faceFactor1_.correctBoundaryConditions();
    faceFactor2_.correctBoundaryConditions();
    faceGrowth1_.correctBoundaryConditions();
    faceGrowth2_.correctBoundaryConditions();
    equivalentGrowth_.correctBoundaryConditions();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::implicitMomentSource
(
    const volScalarMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::crystallizationSource
(
    const labelList& momentOrder,
    const label celli
)
{
    const label i = momentOrder[0];
    const label j = momentOrder[1];
    const label k = momentOrder[2];

    scalar source = 0.0;
    if (nucleation_)
    {
        source =
            nucleationModel_->nucleationSource(i, celli)
           *std::pow(z1Nucleation_.value(), j)
           *std::pow(z2Nucleation_.value(), k);
    }

    if (!growth_)
    {
        return source;
    }

    scalar p1, p2, G1, G2, GL, q1, q2;
    kineticState(celli, p1, p2, G1, G2, GL, q1, q2);
    (void)p1;
    (void)p2;
    (void)G1;
    (void)G2;

    const volScalarMomentFieldSet& moments = quadrature_.moments();

    if (i > 0 && GL > 0)
    {
        labelList lower(momentOrder);
        lower[0] -= 1;
        source += i*GL*max(moments(lower)[celli], scalar(0));
    }

    const scalar incorporationCoefficient =
        3.0*rhop_.value()*shapeFactor_.value()/mRef_.value();

    if (j > 0 && q1 > 0 && GL > 0)
    {
        labelList shifted(momentOrder);
        shifted[0] += 2;
        shifted[1] -= 1;
        source +=
            j*incorporationCoefficient*q1*GL
           *max(moments(shifted)[celli], scalar(0));
    }

    if (k > 0 && q2 > 0 && GL > 0)
    {
        labelList shifted(momentOrder);
        shifted[0] += 2;
        shifted[2] -= 1;
        source +=
            k*incorporationCoefficient*q2*GL
           *max(moments(shifted)[celli], scalar(0));
    }

    return source;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::updateCellMomentSource(const label celli)
{
    crystallizationSourceScale_[celli] = 1.0;
    if (!speciesCoupled_ || (!nucleation_ && !growth_))
    {
        return;
    }

    const scalar invDt =
        1.0/max(phi_.mesh().time().deltaTValue(), scalar(SMALL));
    const scalar hostRate =
        rhop_.value()*shapeFactor_.value()
       *max
        (
            crystallizationSource(labelList({3, 0, 0}), celli),
            scalar(0)
        );
    const scalar impurity1Rate =
        mRef_.value()
       *max
        (
            crystallizationSource(labelList({0, 1, 0}), celli),
            scalar(0)
        );
    const scalar impurity2Rate =
        mRef_.value()
       *max
        (
            crystallizationSource(labelList({0, 0, 1}), celli),
            scalar(0)
        );

    scalar scale = 1.0;
    if (hostRate > SMALL)
    {
        scale = min(scale, 0.9*max(C_[celli], scalar(0))*invDt/hostRate);
    }
    if (impurity1Rate > SMALL)
    {
        scale =
            min(scale, 0.9*max(Ci1_[celli], scalar(0))*invDt/impurity1Rate);
    }
    if (impurity2Rate > SMALL)
    {
        scale =
            min(scale, 0.9*max(Ci2_[celli], scalar(0))*invDt/impurity2Rate);
    }
    crystallizationSourceScale_[celli] = max(scale, scalar(0));
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{
    (void)quadrature;
    (void)environment;
    return
        crystallizationSourceScale_[celli]
       *crystallizationSource(momentOrder, celli);
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::solveMomentSources() const
{
    return (nucleation_ || growth_) && odeType::solveSources_;
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::solveMomentOde() const
{
    return odeType::solveOde_;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::updateStatistics()
{
    const volScalarMomentFieldSet& moments = quadrature_.moments();
    const volScalarMoment& m000 = moments(labelList({0, 0, 0}));
    const volScalarMoment& m100 = moments(labelList({1, 0, 0}));
    const volScalarMoment& m200 = moments(labelList({2, 0, 0}));
    const volScalarMoment& m300 = moments(labelList({3, 0, 0}));
    const volScalarMoment& m010 = moments(labelList({0, 1, 0}));
    const volScalarMoment& m001 = moments(labelList({0, 0, 1}));

    forAll(L10_, celli)
    {
        L10_[celli] = m100[celli]/max(m000[celli], scalar(SMALL));
        L32_[celli] = m300[celli]/max(m200[celli], scalar(SMALL));
        meanZ1_[celli] = m010[celli]/max(m000[celli], scalar(SMALL));
        meanZ2_[celli] = m001[celli]/max(m000[celli], scalar(SMALL));

        const scalar hostMass =
            rhop_.value()*shapeFactor_.value()*max(m300[celli], scalar(0));
        const scalar impurity1Mass =
            mRef_.value()*max(m010[celli], scalar(0));
        const scalar impurity2Mass =
            mRef_.value()*max(m001[celli], scalar(0));
        const scalar totalMass = hostMass + impurity1Mass + impurity2Mass;

        solidCGM1Ppm_[celli] =
            1.0e6*impurity1Mass/max(totalMass, scalar(SMALL));
        solidCGM2Ppm_[celli] =
            1.0e6*impurity2Mass/max(totalMass, scalar(SMALL));
        totalImpurityPpm_[celli] =
            solidCGM1Ppm_[celli] + solidCGM2Ppm_[celli];
        crystalPurityPercent_[celli] =
            100.0*hostMass/max(totalMass, scalar(SMALL));
    }

    L10_.correctBoundaryConditions();
    L32_.correctBoundaryConditions();
    meanZ1_.correctBoundaryConditions();
    meanZ2_.correctBoundaryConditions();
    solidCGM1Ppm_.correctBoundaryConditions();
    solidCGM2Ppm_.correctBoundaryConditions();
    totalImpurityPpm_.correctBoundaryConditions();
    crystalPurityPercent_.correctBoundaryConditions();
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::calcSpeciesTransfer
(
    const scalarField& m300Before,
    const scalarField& m010Before,
    const scalarField& m001Before
)
{
    SYact_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);
    SI1act_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);
    SI2act_ = dimensionedScalar("zero", dimDensity/dimTime, 0.0);

    if (!speciesCoupled_ || !solveMomentSources())
    {
        SYact_.correctBoundaryConditions();
        SI1act_.correctBoundaryConditions();
        SI2act_.correctBoundaryConditions();
        return;
    }

    const scalar invDt =
        1.0/max(phi_.mesh().time().deltaTValue(), scalar(SMALL));
    const volScalarMoment& m300 =
        quadrature_.moments()(labelList({3, 0, 0}));
    const volScalarMoment& m010 =
        quadrature_.moments()(labelList({0, 1, 0}));
    const volScalarMoment& m001 =
        quadrature_.moments()(labelList({0, 0, 1}));

    forAll(SYact_, celli)
    {
        const scalar deltaM300 = m300[celli] - m300Before[celli];
        const scalar deltaM010 = m010[celli] - m010Before[celli];
        const scalar deltaM001 = m001[celli] - m001Before[celli];
        SYact_[celli] =
            rhop_.value()*shapeFactor_.value()
           *max(deltaM300*invDt, scalar(0));
        SI1act_[celli] = mRef_.value()*max(deltaM010*invDt, scalar(0));
        SI2act_[celli] = mRef_.value()*max(deltaM001*invDt, scalar(0));
    }

    SYact_.correctBoundaryConditions();
    SI1act_.correctBoundaryConditions();
    SI2act_.correctBoundaryConditions();
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::realizableCo() const
{
    return momentAdvection_->realizableCo();
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::CoNum() const
{
    return 0.0;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::solve()
{
    crystallizationSourceScale_.setSize(phi_.mesh().nCells(), 1.0);
    crystallizationSourceScale_ = 1.0;

    momentAdvection_->update();
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

    quadrature_.updateQuadrature();
    const scalarField m300Before
    (
        quadrature_.moments()(labelList({3, 0, 0})).primitiveField()
    );
    const scalarField m010Before
    (
        quadrature_.moments()(labelList({0, 1, 0})).primitiveField()
    );
    const scalarField m001Before
    (
        quadrature_.moments()(labelList({0, 0, 1})).primitiveField()
    );

    if (solveMomentSources())
    {
        explicitMomentSource();
    }

    updateKineticDiagnostics();
    updateStatistics();
    calcSpeciesTransfer(m300Before, m010Before, m001Before);
}


bool Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_.subDict(type() + "Coeffs")
    );
    return true;
}


void Foam::PDFTransportModels::populationBalanceModels::
crysDualImpurityPopulationBalance::correctKinetics()
{
    updateKineticDiagnostics();
}

// ************************************************************************* //

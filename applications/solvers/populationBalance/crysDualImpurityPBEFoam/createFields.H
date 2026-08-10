#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

Switch solveFlow(pimple.dict().getOrDefault<bool>("solveFlow", true));
const word temperatureMode
(
    pimple.dict().getOrDefault<word>("temperatureMode", "transport")
);

autoPtr<Function1<scalar>> prescribedTemperature;
if (temperatureMode == "prescribed")
{
    prescribedTemperature =
        Function1<scalar>::New("prescribedTemperature", pimple.dict(), &runTime);
}
else if (temperatureMode != "transport")
{
    FatalIOErrorInFunction(pimple.dict())
        << "Unknown temperatureMode '" << temperatureMode
        << "'. Valid options are transport and prescribed."
        << exit(FatalIOError);
}

singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->nu()
);

#include "createMRF.H"
#include "createFvOptions.H"

Info<< "Reading temperature, solute and two impurity fields\n" << endl;
volScalarField T
(
    IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
volScalarField C
(
    IOobject("C", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
volScalarField Ci1
(
    IOobject("Ci1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
volScalarField Ci2
(
    IOobject("Ci2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word populationBalanceModelName
(
    populationBalanceProperties.lookup("populationBalanceModel")
);
if (populationBalanceModelName != "crysDualImpurityPopulationBalance")
{
    FatalIOErrorInFunction(populationBalanceProperties)
        << "crysDualImpurityPBEFoam requires populationBalanceModel "
        << "crysDualImpurityPopulationBalance, but found "
        << populationBalanceModelName << exit(FatalIOError);
}

const dictionary& crystalCoeffs =
    populationBalanceProperties.subDict(populationBalanceModelName + "Coeffs");

dimensionedScalar Pr
(
    "Pr", dimless, crystalCoeffs.getOrDefault<scalar>("Pr", 1.0)
);
dimensionedScalar Prt
(
    "Prt", dimless, crystalCoeffs.getOrDefault<scalar>("Prt", 1.0)
);
dimensionedScalar molarDiff
(
    "molarDiff",
    dimViscosity,
    crystalCoeffs.getOrDefault<scalar>("molarDiff", 1.0e-9)
);
dimensionedScalar Sct
(
    "Sct", dimless, crystalCoeffs.getOrDefault<scalar>("Sct", 0.7)
);

autoPtr<solutionSaturationModel> saturationModel
(
    solutionSaturationModel::New
    (
        crystalCoeffs.subDict("saturationModel"),
        mesh
    )
);

volScalarField Csat(saturationModel->Csat(T));
volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (C - Csat)/max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL)),
    fvPatchFieldBase::zeroGradientType()
);
sigma.correctBoundaryConditions();

uniformDimensionedScalarField rhop
(
    IOobject
    (
        "rhop",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dimensionedScalar::getOrDefault
    (
        "rhop",
        crystalCoeffs,
        dimDensity,
        2338.0
    )
);
rhop.store();

autoPtr
<
    PDFTransportModels::populationBalanceModels::
    crysDualImpurityPopulationBalance
> populationBalance
(
    new PDFTransportModels::populationBalanceModels::
        crysDualImpurityPopulationBalance
    (
        "populationBalance",
        crystalCoeffs,
        phi
    )
);

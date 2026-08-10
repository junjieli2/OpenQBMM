{
    volScalarField::Internal Su
    (
        IOobject("Su.twoPopulation", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0.0)
    );

    if (populationBalance->speciesCoupled())
    {
        const volScalarField::Internal& SYact =
            populationBalance->SYact();

        forAll(Su, celli)
        {
            Su[celli] = -SYact[celli];
        }
    }

    volScalarField speciesDiffusivity
    (
        IOobject
        (
            "speciesDiffusivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        turbulence->nut()/Sct + molarDiff
    );
    speciesDiffusivity.correctBoundaryConditions();

    fv::EulerDdtScheme<scalar> speciesDdt(mesh);

    fvScalarMatrix CEqn
    (
        speciesDdt.fvmDdt(C)
      + fvm::div(phi, C)
      - fvm::laplacian(speciesDiffusivity, C)
     ==
        Su
      + fvOptions(C)
    );

    fvOptions.constrain(CEqn);
    CEqn.solve();
    fvOptions.correct(C);
    C.correctBoundaryConditions();

    forAll(C, celli)
    {
        if
        (
            !std::isfinite(C[celli])
         || C[celli] < -ROOTSMALL*max(mag(C[celli]), scalar(1))
        )
        {
            FatalErrorInFunction
                << "Invalid solute concentration C=" << C[celli]
                << " in cell " << celli << exit(FatalError);
        }
    }

    C.clamp_min(dimensionedScalar("zeroC", dimDensity, 0.0));
    C.correctBoundaryConditions();

    Csat = saturationModel->Csat(T);
    Csat.correctBoundaryConditions();
    sigma =
        (C - Csat)
       /max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL));
    sigma.correctBoundaryConditions();

    Info<< "C equation solved. Min/max C = "
        << gMin(C) << " / " << gMax(C) << endl;
}

{
    volScalarField::Internal Su
    (
        IOobject("Su_crys", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
    );

    if (populationBalance->speciesCoupled())
    {
        const volScalarField::Internal& SYact = populationBalance->SYact();

        forAll(Su, celli)
        {
            // SYact is the exact solid-phase gain inferred from delta(mu30).
            // Keeping it explicit preserves the equal-and-opposite liquid
            // transfer even when advection and diffusion are present.
            Su[celli] = -SYact[celli];
        }
    }

    // Solute and impurity use the same molecular and turbulent diffusivity.
    // Keep this local so alphat retains its thermal-diffusivity meaning.
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

    // Force first-order Euler for dissolved species. SYact is reconstructed
    // from the solid moment increment over exactly this time interval, so the
    // same Euler interval is required for equal-and-opposite phase transfer.
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

    forAll(C.primitiveField(), celli)
    {
        const scalar value = C[celli];
        const scalar localScale = max(mag(value), scalar(1));

        if (!std::isfinite(value) || value < -ROOTSMALL*localScale)
        {
            FatalErrorInFunction
                << "The coupled solute update produced invalid C=" << value
                << " in cell " << celli << ". Reduce deltaT or use a more "
                << "bounded transport discretization."
                << exit(FatalError);
        }
    }

    Csat = saturationModel->Csat(T);
    Csat.correctBoundaryConditions();
    dimensionedScalar minC("minC", dimDensity, 0.0);
    dimensionedScalar maxC("maxC", dimDensity, GREAT);
    C.clamp_range(minC, maxC);
    C.correctBoundaryConditions();

    sigma =
        (C - Csat)
       /max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL));
    sigma.correctBoundaryConditions();

    Info<< "C equation solved. Min/max C = " << gMin(C) << " / " << gMax(C)
        << endl;
}

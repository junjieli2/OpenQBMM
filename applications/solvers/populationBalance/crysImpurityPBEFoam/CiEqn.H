{
    volScalarField::Internal SuImpurity
    (
        IOobject("Su_impurity", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
    );

    if (populationBalance->speciesCoupled())
    {
        const volScalarField::Internal& SIact = populationBalance->SIact();

        forAll(SuImpurity, celli)
        {
            // SIact is the exact solid-phase gain inferred from delta(mu01).
            SuImpurity[celli] = -SIact[celli];
        }
    }

    // Re-evaluate in case turbulence was corrected; coefficients match C.
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

    // Use the same fixed Euler interval as the solid-to-liquid transfer.
    fv::EulerDdtScheme<scalar> speciesDdt(mesh);

    fvScalarMatrix CiEqn
    (
        speciesDdt.fvmDdt(Ci)
      + fvm::div(phi, Ci)
      - fvm::laplacian(speciesDiffusivity, Ci)
     ==
        SuImpurity
      + fvOptions(Ci)
    );

    fvOptions.constrain(CiEqn);
    CiEqn.solve();
    fvOptions.correct(Ci);
    Ci.correctBoundaryConditions();

    forAll(Ci.primitiveField(), celli)
    {
        const scalar value = Ci[celli];
        const scalar localScale = max(mag(value), scalar(1));

        if (!std::isfinite(value) || value < -ROOTSMALL*localScale)
        {
            FatalErrorInFunction
                << "The coupled impurity update produced invalid Ci="
                << value << " in cell " << celli << ". Reduce deltaT or use "
                << "a more bounded transport discretization."
                << exit(FatalError);
        }
    }

    dimensionedScalar minCi("minCi", dimDensity, 0.0);
    dimensionedScalar maxCi("maxCi", dimDensity, GREAT);
    Ci.clamp_range(minCi, maxCi);
    Ci.correctBoundaryConditions();

    Info<< "Ci equation solved. Min/max Ci = " << gMin(Ci)
        << " / " << gMax(Ci) << endl;
}

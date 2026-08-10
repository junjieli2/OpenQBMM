{
    const volScalarField::Internal& SI1act = populationBalance->SI1act();
    const volScalarField::Internal& SI2act = populationBalance->SI2act();

    auto solveImpurity =
    [&]
    (
        volScalarField& impurity,
        const volScalarField::Internal& solidGain
    )
    {
        volScalarField::Internal source
        (
            IOobject
            (
                "Su_" + impurity.name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
        );

        if (populationBalance->speciesCoupled())
        {
            forAll(source, celli)
            {
                source[celli] = -solidGain[celli];
            }
        }

        volScalarField speciesDiffusivity
        (
            IOobject
            (
                "speciesDiffusivity." + impurity.name(),
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

        fvScalarMatrix impurityEqn
        (
            speciesDdt.fvmDdt(impurity)
          + fvm::div(phi, impurity)
          - fvm::laplacian(speciesDiffusivity, impurity)
         == source + fvOptions(impurity)
        );

        fvOptions.constrain(impurityEqn);
        impurityEqn.solve();
        fvOptions.correct(impurity);
        impurity.correctBoundaryConditions();

        forAll(impurity.primitiveField(), celli)
        {
            const scalar value = impurity[celli];
            const scalar localScale = max(mag(value), scalar(1));
            if (!std::isfinite(value) || value < -ROOTSMALL*localScale)
            {
                FatalErrorInFunction
                    << "The coupled impurity update produced invalid "
                    << impurity.name() << '=' << value << " in cell "
                    << celli << ". Reduce deltaT or use bounded transport."
                    << exit(FatalError);
            }
        }

        dimensionedScalar minimum("minimum", dimDensity, 0.0);
        dimensionedScalar maximum("maximum", dimDensity, GREAT);
        impurity.clamp_range(minimum, maximum);
        impurity.correctBoundaryConditions();

        Info<< impurity.name() << " equation solved. Min/max = "
            << gMin(impurity) << " / " << gMax(impurity) << endl;
    };

    solveImpurity(Ci1, SI1act);
    solveImpurity(Ci2, SI2act);
}

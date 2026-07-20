{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);
    T.correctBoundaryConditions();

    // Keep the PBE driving force synchronized with the just-updated
    // temperature before evaluating crystal source terms.
    Csat = saturationModel->Csat(T);
    Csat.correctBoundaryConditions();
    sigma =
        (C - Csat)
       /max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL));
    sigma.correctBoundaryConditions();
}

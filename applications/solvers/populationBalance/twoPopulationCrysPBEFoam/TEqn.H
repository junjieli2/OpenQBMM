{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);
    fv::EulerDdtScheme<scalar> thermalDdt(mesh);

    fvScalarMatrix TEqn
    (
        thermalDdt.fvmDdt(T)
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
}

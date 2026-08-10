{
    // Keep the PBE driving force synchronized with the current temperature
    // before evaluating crystal source terms.
    Csat = saturationModel->Csat(T);
    Csat.correctBoundaryConditions();
    sigma =
        (C - Csat)
       /max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL));
    sigma.correctBoundaryConditions();
}

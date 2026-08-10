volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA +=
        MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
}

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

constrainPressure(p, U, phiHbyA, rAtU(), MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

fvc::correctUf(Uf, U, phi);
fvc::makeRelative(phi, U);

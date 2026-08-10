/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenQBMM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    crysDualImpurityPBEFoam

Description
    Transient incompressible flow solver coupled to the reduced
    (L,z1,z2) crystal population balance and two dissolved impurities.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "EulerDdtScheme.H"
#include "fvcSmooth.H"
#include "crysDualImpurityPopulationBalance.H"
#include "solutionSaturationModel.H"
#include "Function1.H"
#include <cmath>

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient incompressible flow with reduced L-z1-z2 crystal PBE"
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (LTS)
    {
        FatalErrorInFunction
            << "crysDualImpurityPBEFoam requires a global time step for "
            << "conservative solid/liquid transfer." << exit(FatalError);
    }

    const word velocityDdtScheme(mesh.ddtScheme("ddt(" + U.name() + ')'));
    if (mesh.dynamic() && velocityDdtScheme != "Euler")
    {
        FatalErrorInFunction
            << "A dynamic mesh requires ddt(U) Euler; selected '"
            << velocityDdtScheme << "'." << exit(FatalError);
    }

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (solveFlow && pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();
                if (mesh.changing())
                {
                    MRF.update();
                    if (correctPhi)
                    {
                        phi = mesh.Sf() & Uf();
                        #include "../crysImpurityPBEFoam/correctPhi.H"
                        fvc::makeRelative(phi, U);
                    }
                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "../crysImpurityPBEFoam/UEqn.H"
            while (pimple.correct())
            {
                #include "../crysImpurityPBEFoam/pEqn.H"
            }
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        if (temperatureMode == "transport")
        {
            #include "../crysImpurityPBEFoam/TEqn.H"
        }
        else
        {
            const scalar imposedTemperature =
                prescribedTemperature->value(runTime.value());
            if (!std::isfinite(imposedTemperature) || imposedTemperature <= 0)
            {
                FatalErrorInFunction
                    << "prescribedTemperature must be finite and positive."
                    << exit(FatalError);
            }
            T = dimensionedScalar
            (
                "prescribedTemperature",
                dimTemperature,
                imposedTemperature
            );
            T.correctBoundaryConditions();
        }

        #include "../crysImpurityPBEFoam/updateSaturation.H"
        populationBalance->solve();
        #include "../crysImpurityPBEFoam/CEqn.H"
        #include "CiEqns.H"
        populationBalance->correctKinetics();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //

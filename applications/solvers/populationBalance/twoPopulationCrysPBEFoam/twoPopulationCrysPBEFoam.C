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
#include "twoPopulationBalanceSystem.H"
#include "solutionSaturationModel.H"
#include "Function1.H"

#include <cmath>

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient incompressible flow, temperature, solute and "
        "two-population crystal PBE solver"
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
            << "twoPopulationCrysPBEFoam requires one global Euler time "
            << "interval; local time stepping is not supported."
            << exit(FatalError);
    }

    const word velocityDdtScheme
    (
        mesh.ddtScheme("ddt(" + U.name() + ')')
    );

    if (mesh.dynamic() && velocityDdtScheme != "Euler")
    {
        FatalErrorInFunction
            << "A dynamic mesh requires ddt(U) Euler so both populations, "
            << "temperature and solute share the same geometry update."
            << exit(FatalError);
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
                        #include "correctPhi.H"
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        if (temperatureMode == "transport")
        {
            #include "TEqn.H"
        }
        else
        {
            const scalar imposedTemperature =
                prescribedTemperature->value(runTime.value());

            if
            (
                !std::isfinite(imposedTemperature)
             || imposedTemperature <= 0
            )
            {
                FatalErrorInFunction
                    << "prescribedTemperature must be finite and positive; "
                    << "received " << imposedTemperature << " K at t="
                    << runTime.value() << exit(FatalError);
            }

            T = dimensionedScalar
            (
                "prescribedTemperature",
                dimTemperature,
                imposedTemperature
            );
            T.correctBoundaryConditions();
        }

        Csat = saturationModel->Csat(T);
        Csat.correctBoundaryConditions();
        sigma =
            (C - Csat)
           /max(Csat, dimensionedScalar("minCsat", dimDensity, SMALL));
        sigma.correctBoundaryConditions();

        nu = turbulence->nu();
        nu.correctBoundaryConditions();

        populationBalance->solve();

        #include "CEqn.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM and is distributed under the GNU General
    Public License, version 3 or later.

Application
    crysImpurityPBEFoam

Description
    Transient incompressible flow solver coupled to a bivariate crystal-size
    and impurity population balance, solute transport and impurity transport.

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
#include "crysImpurityPopulationBalance.H"
#include "solutionSaturationModel.H"

#include <cmath>

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient incompressible flow with bivariate crystal and impurity "
        "population-balance transport"
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
            << "crysImpurityPBEFoam requires a single global time step so "
            << "that crystal and dissolved-species transfers remain "
            << "conservative; local time stepping is not supported."
            << exit(FatalError);
    }

    const word velocityDdtScheme
    (
        mesh.ddtScheme("ddt(" + U.name() + ')')
    );

    if (mesh.dynamic() && velocityDdtScheme != "Euler")
    {
        FatalErrorInFunction
            << "A dynamic mesh requires ddt(U) Euler in "
            << "system/fvSchemes. The relative flux meshPhi must use the "
            << "same Euler geometry update as the bivariate moments, C and "
            << "Ci, but the selected ddt(U) scheme is '"
            << velocityDdtScheme << "'."
            << exit(FatalError);
    }

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

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

        #include "TEqn.H"

        populationBalance->solve();

        #include "CEqn.H"
        #include "CiEqn.H"

        populationBalance->correctAdsorption();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

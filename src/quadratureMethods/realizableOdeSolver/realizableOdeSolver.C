/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "realizableOdeSolver.H"
#include "PstreamReduceOps.H"

#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
Foam::realizableOdeSolver<momentType, nodeType>::realizableOdeSolver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    ATol_(readScalar(dict.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(dict.subDict("odeCoeffs").lookup("RTol"))),
    fac_(readScalar(dict.subDict("odeCoeffs").lookup("fac"))),
    facMin_(readScalar(dict.subDict("odeCoeffs").lookup("facMin"))),
    facMax_(readScalar(dict.subDict("odeCoeffs").lookup("facMax"))),
    minLocalDt_(readScalar(dict.subDict("odeCoeffs").lookup("minLocalDt"))),
    localDt_
    (
        IOobject
        (
            "realizableOde:localDt",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        mesh.time().deltaT()
    ),
    localDtAdjustments_(0),
    scaleMoments_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("momentScaling", false)
    ),
    ATolNorm_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("ATolNorm", 1.0e-12)
    ),
    momentScales_(),
    diagOde_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("diagOde", false)
    ),
    diagInterval_
    (
        dict.subDict("odeCoeffs").lookupOrDefault<label>("diagInterval", 1)
    ),
    diagStep_(0),
    nSubStepsMax_(0),
    nSubStepsSum_(0),
    nRejected_(0),
    localDtMin_(0),
    localDtMax_(0),
    errorMax_(0),
    solveSources_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("solveSources", true)
    ),
    solveOde_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("solveOde", true)
    )
{
    // A restarted run must not silently inherit a corrupt or foreign local
    // step field: the cached step controls the ODE sub-stepping and therefore
    // the numerical trajectory.
    if (localDt_.size() != mesh_.nCells())
    {
        FatalIOErrorInFunction(dict)
            << "Field " << localDt_.name() << " has "
            << localDt_.size() << " entries but the mesh has "
            << mesh_.nCells() << " cells." << exit(FatalIOError);
    }

    // Values above the current global step are legal: solve() clamps the
    // cached step with min(localDt_[celli], globalDt). Only non-finite or
    // non-positive entries are unusable.
    label nBad = 0;
    scalar firstBad = 0;
    label firstBadCell = -1;

    forAll(localDt_, celli)
    {
        const scalar value = localDt_[celli];

        if (!std::isfinite(value) || value <= 0.0)
        {
            if (nBad == 0)
            {
                firstBad = value;
                firstBadCell = celli;
            }
            nBad++;
        }
    }

    if (nBad)
    {
        FatalIOErrorInFunction(dict)
            << "Field " << localDt_.name() << " has " << nBad
            << " invalid entries; first at cell " << firstBadCell
            << " with value " << firstBad
            << ". Every entry must be finite and positive."
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
Foam::realizableOdeSolver<momentType, nodeType>::~realizableOdeSolver()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
void Foam::realizableOdeSolver<momentType, nodeType>::solve
(
    quadratureType& quadrature,
    const label enviroment
)
{
    if (!solveSources_)
    {
        return;
    }

    momentFieldSetType& moments(quadrature.moments());
    label nMoments = quadrature.nMoments();
    scalar globalDt = mesh_.time().deltaT().value();
    const labelListList& momentOrders = quadrature.momentOrders();

    // Reset the local solver diagnostics recorded by this call
    nSubStepsMax_ = 0;
    nSubStepsSum_ = 0;
    nRejected_ = 0;
    localDtMin_ = GREAT;
    localDtMax_ = 0.0;
    errorMax_ = 0.0;

    if (scaleMoments_)
    {
        initialiseMomentScales(moments, nMoments);
    }

    //- Use Euler explicit to update moments due to sources
    if (!solveOde_)
    {
        forAll(moments[0], celli)
        {
            updateCellMomentSource(celli);

            scalarList momentSources(nMoments, Zero);

            forAll(moments, mi)
            {
                const labelList& order = momentOrders[mi];

                momentSources[mi] = cellMomentSource
                (
                    order,
                    celli,
                    quadrature,
                    enviroment
                );
            }

            // Apply all components together so every source evaluation sees
            // the same frozen moment state.
            forAll(moments, mi)
            {
                moments[mi][celli] += globalDt*momentSources[mi];
            }

            quadrature.updateLocalQuadrature(celli, true);
            quadrature.updateLocalMoments(celli);
        }

        forAll(moments, mi)
        {
            moments[mi].correctBoundaryConditions();
        }

        quadrature.updateBoundaryQuadrature();

        return;
    }

    Info << "Solving source terms in realizable ODE solver." << endl;

    forAll(moments[0], celli)
    {
        //Info << "OLD MOMENTS" << moments << endl;

        // Moments are the primary ODE state. Invert them to obtain nodes for
        // source evaluation, but do not overwrite them with moments rebuilt
        // from an ill-conditioned conditional quadrature: that projection can
        // destroy exact source invariants such as paired host/impurity mass.
        // Store the primary moments to recover from a failed step.
        quadrature.updateLocalQuadrature(celli, true, false);

        scalarList oldMoments(nMoments, Zero);

        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }

        //Info << "Old moments: " << oldMoments << endl;

        //- Local time
        scalar localT(0);

        // Initialize the local step
        scalar localDt = min(localDt_[celli], globalDt);

        localDtMin_ = min(localDtMin_, localDt);
        localDtMax_ = max(localDtMax_, localDt);

        label cellSubSteps = 0;

        // Initialize RK parameters
        scalarList k1(nMoments, Zero);
        scalarList k2(nMoments, Zero);
        scalarList k3(nMoments, Zero);

        // Flag to indicate if the time step is complete
        bool timeComplete = false;

        // Check realizability of intermediate moment sets
        bool realizableUpdate1 = false;
        bool realizableUpdate2 = false;
        bool realizableUpdate3 = false;

        scalarList diff23(nMoments, Zero);
        label nItt = 0;

        while (!timeComplete)
        {
            do
            {
                nItt++;
                cellSubSteps++;

                // First intermediate update
                bool nullSource =  true;

                updateCellMomentSource(celli);

                forAll(k1, mi)
                {
                    const labelList& order = momentOrders[mi];

                    k1[mi] =
                        localDt*cellMomentSource
                        (
                            order,
                            celli,
                            quadrature,
                            enviroment
                        );

                    if (mag(k1[mi]) > SMALL)
                    {
                        nullSource = false;
                    }
                }

                forAll(moments, mi)
                {
                    moments[mi][celli] = oldMoments[mi] + k1[mi];
                }

                realizableUpdate1 =
                        quadrature.updateLocalQuadrature(celli, false, false);


                if (nullSource)
                {
                    // With a vanishing source every stage update is zero and so
                    // is the embedded error. diff23 persists across cells, so it
                    // must be cleared explicitly; otherwise the error test below
                    // would score this cell with another cell's estimate.
                    forAll(diff23, mi)
                    {
                        diff23[mi] = Zero;
                    }

                    break;
                }

                // Second moment update
                updateCellMomentSource(celli);

                forAll(k2, mi)
                {
                    const labelList& order = momentOrders[mi];

                    k2[mi] =
                        localDt*cellMomentSource
                        (
                            order,
                            celli,
                            quadrature,
                            enviroment
                        );

                }

                forAll(moments, mi)
                {
                    moments[mi][celli] =
                        oldMoments[mi] + (k1[mi] + k2[mi])/4.0;
                }

                realizableUpdate2 =
                    quadrature.updateLocalQuadrature(celli, false, false);


                // Third moment update
                updateCellMomentSource(celli);

                forAll(k3, mi)
                {
                    const labelList& order = momentOrders[mi];

                    k3[mi] =
                        localDt*cellMomentSource
                        (
                            order,
                            celli,
                            quadrature,
                            enviroment
                        );

                    diff23[mi] = (2.0*k3[mi] - k1[mi] - k2[mi])/3.0;
                }


                forAll(moments, mi)
                {
                    moments[mi][celli] =
                        oldMoments[mi]
                      + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;
                }

                realizableUpdate3 =
                    quadrature.updateLocalQuadrature(celli, false, false);


                if
                (
                    realizableUpdate1
                 && realizableUpdate2
                 && realizableUpdate3
                 && !acceptMomentUpdate(celli)
                )
                {
                    // The trial step is realizable but the model rejected it,
                    // for example because it would consume more solute than the
                    // current split step has available. Route it through the
                    // same recovery as a realizability failure.
                    realizableUpdate3 = false;
                }

                if
                (
                    !realizableUpdate1
                 || !realizableUpdate2
                 || !realizableUpdate3
                )
                {
                    // Avoid spamming the terminal on repeated rejections
                    if (localDtAdjustments_ == 0)
                    {
                        Info << "Local ODE step rejected, adjusting local "
                             << "timestep." << nl
                             << "This may take a while." << endl;
                    }

                    localDtAdjustments_++;
                    nRejected_++;

                    forAll(oldMoments, mi)
                    {
                        moments[mi][celli] = oldMoments[mi];
                    }

                    // Updating local quadrature with old moments
                    quadrature.updateLocalQuadrature(celli, true, false);

                    localDt /= 2.0;

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minimum local step in realizable ODE"
                            << nl
                            << "    solver. Cannot ensure realizability."
                            << nl
                            << abort(FatalError);
                    }
                }
            }
            while
            (
                !realizableUpdate1
             || !realizableUpdate2
             || !realizableUpdate3
            );

            // Initialize error and change
            scalar error(0);
            scalar maxChange(0);

            const bool useScales =
                scaleMoments_ && momentScales_.size() == nMoments;

            for (label mi = 0; mi < nMoments; mi++)
            {
                // Calculate the scaling factor. The relative term is always
                // present; the absolute term is either the raw ATol_ or, with
                // momentScaling enabled, ATolNorm_ referred to the magnitude of
                // that moment. A single raw ATol cannot serve a moment set
                // spanning many orders of magnitude: for the small moments it
                // exceeds every admissible step error, which silently disables
                // the embedded error control.
                scalar scalei =
                    max
                    (
                        mag(moments[mi][celli]), mag(oldMoments[mi])
                    )*RTol_
                  + (useScales ? momentScales_[mi]*ATolNorm_ : ATol_);

                scalei = max(scalei, VSMALL);

                // Update the error
                error += sqr(diff23[mi]/scalei);

                // Update the maximum change in moments
                maxChange
                    = max(maxChange, mag(moments[mi][celli] - oldMoments[mi]));
            }

            error = sqrt(error/nMoments);
            errorMax_ = max(errorMax_, error);

            // Fac_/cbrt(error) is the step-size factor. Flooring the radicand
            // keeps a vanishing or non-finite estimate from trapping a
            // floating-point divide instead of taking the normal accept path.
            const scalar errRoot = pow(max(error, VSMALL), 1.0/3.0);

            if (error < SMALL)
            {
                // An embedded error of zero means that this local step was
                // accepted, not that the complete global time interval was
                // covered. This distinction matters when deltaT increases
                // after a smaller local step was cached (for example, for a
                // constant nucleation source).
                localT += localDt;

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                const scalar remaining =
                    max(globalDt - localT, scalar(0));

                if (remaining <= SMALL*max(globalDt, scalar(1)))
                {
                    timeComplete = true;
                    localT = Zero;
                    break;
                }

                // With vanishing estimated error, try the complete remaining
                // interval. A nonlinear source will still be rejected and
                // reduced by the normal error-control branch if necessary.
                localDt = remaining;
                localDt_[celli] = localDt;
            }
            else if (maxChange < SMALL)
            {
                WarningInFunction
                    << "The maximum change in moments is small, "
                    << "but error is not.\n"
                    << nl
                    << "Error: " << error << nl
                    << "Max. change: " << maxChange << nl
                    << nl
                    << "\nThis may indicate a problem with the "
                    << "realizable ODE solver." << endl;

                // The current change is negligible in absolute terms, but it
                // still only covers localDt. Accept it and try the remaining
                // global interval so a previously cached small step cannot
                // silently truncate the physical time advance.
                localT += localDt;

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                const scalar remaining =
                    max(globalDt - localT, scalar(0));

                if (remaining <= SMALL*max(globalDt, scalar(1)))
                {
                    timeComplete = true;
                    localT = Zero;
                    break;
                }

                localDt = remaining;
                localDt_[celli] = localDt;
            }
            else if (error < 1)
            {
                localT += localDt;
                localDt *= min(facMax_, max(facMin_, fac_/errRoot));
                scalar maxLocalDt = max(globalDt - localT, scalar(0));
                localDt = min(maxLocalDt, localDt);

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                if (localDt == 0.0)
                {
                    timeComplete = true;
                    localT = Zero;
                    break;
                }

                localDt_[celli] = localDt;
            }
            else
            {
                localDt *=
                    min(scalar(1), max(facMin_, fac_/errRoot));

                nRejected_++;

                forAll(oldMoments, mi)
                {
                    moments[mi][celli] = oldMoments[mi];
                }

                // Updating local quadrature with old moments
                quadrature.updateLocalQuadrature(celli, true, false);
            }
        }

        nSubStepsSum_ += cellSubSteps;
        nSubStepsMax_ = max(nSubStepsMax_, cellSubSteps);
    }

    forAll(moments, mi)
    {
        moments[mi].correctBoundaryConditions();
    }

    quadrature.updateBoundaryQuadrature();

    diagStep_++;

    if (diagOde_ && diagInterval_ > 0 && (diagStep_ % diagInterval_) == 0)
    {
        label nCells = moments[0].size();
        reduce(nCells, sumOp<label>());
        reduce(nSubStepsMax_, maxOp<label>());
        reduce(nSubStepsSum_, sumOp<label>());
        reduce(nRejected_, sumOp<label>());
        reduce(localDtMin_, minOp<scalar>());
        reduce(localDtMax_, maxOp<scalar>());
        reduce(errorMax_, maxOp<scalar>());

        Info<< "realizableOde: t = " << mesh_.time().timeName()
            << " globalDt = " << globalDt
            << " maxSubSteps = " << nSubStepsMax_
            << " sumSubSteps = " << nSubStepsSum_
            << " nCells = " << nCells
            << " rejected = " << nRejected_
            << " localDtMin = " << localDtMin_
            << " localDtMax = " << localDtMax_
            << " errorMax = " << errorMax_ << endl;
    }
}


template<class momentType, class nodeType>
bool Foam::realizableOdeSolver<momentType, nodeType>::acceptMomentUpdate
(
    const label celli
)
{
    return true;
}


template<class momentType, class nodeType>
void Foam::realizableOdeSolver<momentType, nodeType>::initialiseMomentScales
(
    const momentFieldSetType& moments,
    const label nMoments
) const
{
    if (momentScales_.size() == nMoments)
    {
        return;
    }

    momentScales_.setSize(nMoments, 1.0);

    // A single fixed scale per moment order: the largest magnitude that moment
    // attains anywhere on the mesh at the time the scales are first needed.
    // Deriving it once keeps the error control from drifting with the solution,
    // which would make the accepted step history-dependent.
    for (label mi = 0; mi < nMoments; mi++)
    {
        scalar cellMax = 0.0;

        forAll(moments[mi], celli)
        {
            const scalar value = mag(moments[mi][celli]);

            if (std::isfinite(value))
            {
                cellMax = max(cellMax, value);
            }
        }

        reduce(cellMax, maxOp<scalar>());

        momentScales_[mi] = max(cellMax, VSMALL);
    }

    Info<< "realizableOde: moment scales initialised to " << momentScales_
        << endl;
}


template<class momentType, class nodeType>
void Foam::realizableOdeSolver<momentType, nodeType>
::read(const dictionary& dict)
{
    const dictionary& odeDict = dict.subDict("odeCoeffs");
    solveSources_ = odeDict.lookupOrDefault<Switch>("solveSources", true);
    solveOde_ = odeDict.lookupOrDefault<Switch>("solveOde", true);
    scaleMoments_ = odeDict.lookupOrDefault<Switch>("momentScaling", false);
    diagOde_ = odeDict.lookupOrDefault<Switch>("diagOde", false);
    diagInterval_ = odeDict.lookupOrDefault<label>("diagInterval", 1);
    ATolNorm_ = odeDict.lookupOrDefault("ATolNorm", 1.0e-12);

    (odeDict.lookup("ATol")) >> ATol_;
    (odeDict.lookup("RTol")) >> RTol_;
    (odeDict.lookup("fac")) >> fac_;
    (odeDict.lookup("facMin")) >> facMin_;
    (odeDict.lookup("facMax")) >> facMax_;
    (odeDict.lookup("minLocalDt")) >> minLocalDt_;
}


// ************************************************************************* //

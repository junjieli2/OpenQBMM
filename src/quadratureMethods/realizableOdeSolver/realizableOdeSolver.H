/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2023 Alberto Passalacqua
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

Class
    Foam::realizableOdeSolver

Description
    Integrates a moment vector in a computational cell using a realizable ODE
    solver.

    The procedure is adapted from the Nguyen et al. (2016) to fit into the 
    OpenFOAM structure. Specifically, the time step is calculated to ensure
    realizability and fit within the global time step of the calling solver.

    Reference:
    \verbatim
        "Solution of population balance equations in applications with 
        fine particles: Mathematical modeling and numerical schemes"
        Nguyen, T.T., Laurent, F., Fox, R.O., Massot, M.
        Issue 325, 2016, pp. 129-156.
    \endverbatim

SourceFiles
    realizableOdeSolver.C

\*---------------------------------------------------------------------------*/

#ifndef realizableOdeSolver_H
#define realizableOdeSolver_H

#include "quadratureApproximation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class realizableOdeSolver Declaration
\*---------------------------------------------------------------------------*/

template<class momentType, class nodeType>
class realizableOdeSolver
{

    typedef momentFieldSet<momentType, nodeType> momentFieldSetType;
    typedef quadratureApproximation<momentType, nodeType> quadratureType;

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Absolute tolerance for realizable ODE solver
        scalar ATol_;

        //- Relative tolerance for realizable ODE solver
        scalar RTol_;

        //- Time-step change factor
        scalar fac_;

        //- Minimum time-step change factor
        scalar facMin_;

        //- Maximum time-step change factor
        scalar facMax_;

        //- Minimum local step in realizable ODE solver
        scalar minLocalDt_;

        //- Stored time step values
        volScalarField localDt_;

        //- Number of times localDt has been adjusted
        label localDtAdjustments_;

        //- Solve the embedded error test in normalised moment variables
        Switch scaleMoments_;

        //- Dimensionless absolute tolerance used by the normalised error test
        scalar ATolNorm_;

        //- Per-moment normalisation scales (empty until first solve())
        mutable scalarList momentScales_;

        //- Report per-step local solver statistics
        Switch diagOde_;

        //- Report every this many global steps when diagOde is on
        label diagInterval_;

        //- Global steps seen since construction, used by diagInterval_
        mutable label diagStep_;

        //- Maximum number of local steps taken by any cell in the last solve
        mutable label nSubStepsMax_;

        //- Total number of local steps taken in the last solve
        mutable label nSubStepsSum_;

        //- Number of rejected local steps in the last solve
        mutable label nRejected_;

        //- Smallest local step attempted in the last solve
        mutable scalar localDtMin_;

        //- Largest local step attempted in the last solve
        mutable scalar localDtMax_;

        //- Largest embedded error estimate seen in the last solve
        mutable scalar errorMax_;


    // Private member functions

        //- Populate momentScales_ from the current moment magnitudes
        void initialiseMomentScales
        (
            const momentFieldSetType& moments,
            const label nMoments
        ) const;


protected:

    // Protected data

        //- Switch to activate moment update due to source terms
        Switch solveSources_;

        //- Switch to activate the realizable ODE solver
        Switch solveOde_;

    // Protected member functions

        //- Calculate source term for moment equation in the specified cell
        virtual void updateCellMomentSource(const label celli) = 0;

        //- Additional acceptance test for a completed local step.
        //  Derived models use this to reject a step that is realizable but
        //  unphysical, for example one that would consume more solute than the
        //  current global step has available. Returning false halves the local
        //  step exactly as a realizability failure does.
        virtual bool acceptMomentUpdate(const label celli);

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const labelList& momentOrder,
            const label celli,
            const quadratureType& nodes,
            const label enviroment
        ) = 0;


public:

    // Constructors

        //- Construct from components
        realizableOdeSolver
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~realizableOdeSolver();


    // Member Functions

        //- Solve univariate PDF transport equation
        void solve(quadratureType& quadrature, const label enviroment);

        //- Read coefficients from dictionary
        void read(const dictionary& dict);


        // Local solver diagnostics from the most recent solve()

            //- Maximum number of local steps taken by any cell
            label nSubStepsMax() const
            {
                return nSubStepsMax_;
            }

            //- Total number of local steps taken over all cells
            label nSubStepsSum() const
            {
                return nSubStepsSum_;
            }

            //- Number of rejected local steps
            label nRejected() const
            {
                return nRejected_;
            }

            //- Smallest local step attempted
            scalar localDtMin() const
            {
                return localDtMin_;
            }

            //- Largest local step attempted
            scalar localDtMax() const
            {
                return localDtMax_;
            }

            //- Largest embedded error estimate
            scalar errorMax() const
            {
                return errorMax_;
            }

            //- Per-moment normalisation scales (unit values when disabled)
            const scalarList& momentScales() const
            {
                return momentScales_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if NoRepository
    #include "realizableOdeSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

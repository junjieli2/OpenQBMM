/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is distributed under the GNU General Public License v3 or later.
\*---------------------------------------------------------------------------*/

#include "visibleBoundaryFlux.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(visibleBoundaryFlux, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        visibleBoundaryFlux,
        dictionary
    );
}
}
}


Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::visibleBoundaryFlux
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    amplitude_("amplitude", inv(dimVolume*dimTime), dict),
    startTime_(dict.lookupOrDefault<scalar>("startTime", 0.0)),
    timeScale_(readScalar(dict.lookup("timeScale"))),
    shapeExponent_(dict.lookupOrDefault<label>("shapeExponent", 2)),
    inletMedian_("inletMedian", dimLength, dict),
    inletSigmaLog_(dict.lookupOrDefault<scalar>("inletSigmaLog", 0.0)),
    initialUnresolvedM3_
    (
        readScalar(dict.lookup("initialUnresolvedM3"))
    ),
    fluxRate_
    (
        IOobject
        (
            "visibleBoundaryFluxRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", inv(dimVolume*dimTime), 0.0)
    ),
    unresolvedM3_
    (
        IOobject
        (
            "unresolvedM3",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialUnresolvedM3",
            dimless,
            initialUnresolvedM3_
        )
    ),
    cumulativeTransferredM3_
    (
        IOobject
        (
            "cumulativeVisibleTransferredM3",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    lastTimeIndex_(mesh.time().timeIndex()),
    cumulativeTransferredM3Value_(0.0)
{
    cumulativeTransferredM3Value_ =
        gMax(cumulativeTransferredM3_.primitiveField());

    if (timeScale_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "timeScale must be positive. Received " << timeScale_
            << exit(FatalIOError);
    }

    if (shapeExponent_ < 0 || shapeExponent_ > 20)
    {
        FatalIOErrorInFunction(dict)
            << "shapeExponent must be an integer in [0, 20]. Received "
            << shapeExponent_ << exit(FatalIOError);
    }

    if (inletMedian_.value() <= SMALL || inletSigmaLog_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "inletMedian must be positive and inletSigmaLog must be "
            << "non-negative." << exit(FatalIOError);
    }

    if (initialUnresolvedM3_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "initialUnresolvedM3 must be non-negative. Received "
            << initialUnresolvedM3_ << exit(FatalIOError);
    }

    const scalar totalTransferredM3 =
        amplitude_.value()
       *timeScale_
       *factorial(shapeExponent_)
       *inletRawMoment(3);

    const scalar tolerance =
        1.0e-10*max(initialUnresolvedM3_, scalar(1));

    if (totalTransferredM3 > initialUnresolvedM3_ + tolerance)
    {
        FatalIOErrorInFunction(dict)
            << "The prescribed visible-boundary pulse transfers M3="
            << totalTransferredM3 << ", but initialUnresolvedM3="
            << initialUnresolvedM3_ << ". Increase the reservoir or reduce "
            << "amplitude/timeScale/inlet size." << exit(FatalIOError);
    }

    Info<< "visibleBoundaryFlux: A=" << amplitude_.value()
        << " 1/(m3 s), startTime=" << startTime_
        << " s, timeScale=" << timeScale_
        << " s, shapeExponent=" << shapeExponent_
        << ", inletMedian=" << inletMedian_.value()
        << " m, inletSigmaLog=" << inletSigmaLog_
        << ", total transferred M3=" << totalTransferredM3
        << ", initial unresolved M3=" << initialUnresolvedM3_
        << endl;
}


Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::~visibleBoundaryFlux()
{}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::factorial(const label n) const
{
    scalar value = 1.0;

    for (label i = 2; i <= n; ++i)
    {
        value *= scalar(i);
    }

    return value;
}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::inletRawMoment(const label order) const
{
    const scalar k = scalar(order);

    return
        std::pow(inletMedian_.value(), k)
       *std::exp(0.5*sqr(k*inletSigmaLog_));
}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::numberFlux() const
{
    const scalar elapsed = mesh_.time().value() - startTime_;

    if (elapsed < 0.0)
    {
        return 0.0;
    }

    const scalar x = elapsed/timeScale_;

    return amplitude_.value()*std::pow(x, shapeExponent_)*std::exp(-x);
}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::cumulativeTransferredNumber() const
{
    const scalar elapsed = mesh_.time().value() - startTime_;

    if (elapsed <= 0.0)
    {
        return 0.0;
    }

    const scalar x = elapsed/timeScale_;
    scalar series = 1.0;
    scalar term = 1.0;

    for (label k = 1; k <= shapeExponent_; ++k)
    {
        term *= x/scalar(k);
        series += term;
    }

    const scalar regularized =
        min(max(1.0 - std::exp(-x)*series, scalar(0)), scalar(1));

    return
        amplitude_.value()
       *timeScale_
       *factorial(shapeExponent_)
       *regularized;
}


void
Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::updateDiagnostics(const label celli) const
{
    const label currentTimeIndex = mesh_.time().timeIndex();

    if (currentTimeIndex != lastTimeIndex_)
    {
        cumulativeTransferredM3Value_ +=
            numberFlux()
           *inletRawMoment(3)
           *mesh_.time().deltaTValue();

        lastTimeIndex_ = currentTimeIndex;

        const scalar tolerance =
            1.0e-10*max(initialUnresolvedM3_, scalar(1));

        if
        (
            cumulativeTransferredM3Value_
          > initialUnresolvedM3_ + tolerance
        )
        {
            FatalErrorInFunction
                << "The discretely integrated visible-boundary transfer "
                << "has exhausted the unresolved M3 reservoir: cumulative="
                << cumulativeTransferredM3Value_
                << ", reservoir=" << initialUnresolvedM3_
                << ". Increase initialUnresolvedM3, reduce the pulse, or "
                << "reduce deltaT." << exit(FatalError);
        }
    }

    fluxRate_[celli] = numberFlux();
    cumulativeTransferredM3_[celli] = cumulativeTransferredM3Value_;
    unresolvedM3_[celli] =
        max(initialUnresolvedM3_ - cumulativeTransferredM3Value_, scalar(0));
}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::visibleBoundaryFlux
::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    updateDiagnostics(celli);

    return fluxRate_[celli]*inletRawMoment(momentOrder);
}


// ************************************************************************* //

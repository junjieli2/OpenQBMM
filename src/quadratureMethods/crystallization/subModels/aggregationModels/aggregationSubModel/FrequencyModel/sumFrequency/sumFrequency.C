#include "sumFrequency.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalCollisionFrequencies
{
    defineTypeNameAndDebug(sumFrequency, 0);

    addToRunTimeSelectionTable
    (
        crystalCollisionFrequency,
        sumFrequency,
        dictionary
    );
}
}
}
}

Foam::populationBalanceSubModels::aggregationKernels::
crystalCollisionFrequencies::sumFrequency::sumFrequency
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    crystalCollisionFrequency(dict, mesh),
    componentDicts_(),
    components_()
{
    if (!dict.found("frequencyModels"))
    {
        FatalIOErrorInFunction(dict)
            << "The sum crystalCollisionFrequency requires a "
            << "frequencyModels dictionary containing at least one "
            << "named component." << exit(FatalIOError);
    }

    const dictionary& modelDicts = dict.subDict("frequencyModels");
    componentDicts_.setSize(modelDicts.size());
    components_.setSize(modelDicts.size());

    label modeli = 0;
    forAllConstIter(dictionary, modelDicts, iter)
    {
        componentDicts_.set(modeli, new dictionary(iter().dict()));

        if
        (
            word(componentDicts_[modeli].lookup("crystalCollisionFrequency"))
         == type()
        )
        {
            FatalIOErrorInFunction(componentDicts_[modeli])
                << "A sum frequency cannot contain itself recursively."
                << exit(FatalIOError);
        }

        components_.set
        (
            modeli,
            crystalCollisionFrequency::New
            (
                componentDicts_[modeli],
                mesh
            ).ptr()
        );
        ++modeli;
    }

    if (components_.empty())
    {
        FatalIOErrorInFunction(modelDicts)
            << "frequencyModels must contain at least one component."
            << exit(FatalIOError);
    }
}


Foam::populationBalanceSubModels::aggregationKernels::
crystalCollisionFrequencies::sumFrequency::~sumFrequency()
{}


void Foam::populationBalanceSubModels::aggregationKernels::
crystalCollisionFrequencies::sumFrequency::update()
{
    forAll(components_, modeli)
    {
        components_[modeli].update();
    }
}


Foam::scalar Foam::populationBalanceSubModels::aggregationKernels::
crystalCollisionFrequencies::sumFrequency::beta
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    scalar total = 0.0;

    forAll(components_, modeli)
    {
        const scalar contribution =
            components_[modeli].beta(d1, d2, Ur, celli);

        if (std::isfinite(contribution) && contribution > 0.0)
        {
            total += contribution;
        }
    }

    return std::isfinite(total) ? total : scalar(0);
}

// ************************************************************************* //

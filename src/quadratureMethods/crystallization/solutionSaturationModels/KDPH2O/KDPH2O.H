/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM and is distributed under the
    GNU General Public License, version 3 or later.

Class
    Foam::solutionSaturationModels::KDPH2O

Description
    Temperature-dependent solubility of anhydrous potassium dihydrogen
    phosphate (KDP) in water.

    The published correlation returns kg KDP per kg water:

        r(T_C) = 9.3027e-5*T_C^2 - 9.7629e-5*T_C + 0.2087

    This model converts it to kg KDP per m3 water using rhoSolvent.

\*---------------------------------------------------------------------------*/

#ifndef solutionSaturationModels_KDPH2O_H
#define solutionSaturationModels_KDPH2O_H

#include "solutionSaturationModel.H"

namespace Foam
{
namespace solutionSaturationModels
{

class KDPH2O
:
    public solutionSaturationModel
{
    dimensionedScalar rhoSolvent_;

    scalar saturationConcentration(const scalar temperatureKelvin) const;

public:

    TypeName("KDPH2O");

    KDPH2O(const dictionary& dict, const objectRegistry& db);

    virtual ~KDPH2O() = default;

    virtual tmp<volScalarField> Csat
    (
        const volScalarField& T
    ) const;
};

} // End namespace solutionSaturationModels
} // End namespace Foam

#endif

// ************************************************************************* //

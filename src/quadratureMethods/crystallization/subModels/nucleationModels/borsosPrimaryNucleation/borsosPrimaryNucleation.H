/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenQBMM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenQBMM and is distributed under the GNU General
    Public License, version 3 or later.

Class
    Foam::populationBalanceSubModels::nucleationModels::
    borsosPrimaryNucleation

Description
    Primary-nucleation expression used by Borsos et al. (2016):

        Bp = kp0 exp(-Ep/(R T)) exp(-ke/sqr(log(C/Csat)))

    Nucleation is zero at and below saturation. The returned source for a
    length moment of order q is Bp*pow(Lnucleation, q).

\*---------------------------------------------------------------------------*/

#ifndef borsosPrimaryNucleation_H
#define borsosPrimaryNucleation_H

#include "nucleationModel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

class borsosPrimaryNucleation
:
    public nucleationModel
{
    dimensionedScalar kp0_;
    dimensionedScalar Ep_;
    dimensionedScalar R_;
    scalar ke_;
    scalar nucleationScale_;
    const volScalarField& T_;
    const volScalarField& C_;
    const volScalarField& Csat_;
    volScalarField rate_;
    dimensionedScalar Lnucleation_;

public:

    TypeName("borsosPrimaryNucleation");

    borsosPrimaryNucleation(const dictionary& dict, const fvMesh& mesh);

    virtual ~borsosPrimaryNucleation();

    virtual scalar nucleationSource
    (
        const label& momentOrder,
        const label celli,
        const label environment = 0
    ) const;
};

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::JohnsonJacksonFrictionalStress

Description

SourceFiles
    JohnsonJacksonFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonJacksonFrictionalStress_H
#define JohnsonJacksonFrictionalStress_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JohnsonJacksonFrictionalStress Declaration
\*---------------------------------------------------------------------------*/

class JohnsonJacksonFrictionalStress
:
    public frictionalStressModel
{

public:

    //- Runtime type information
    TypeName("JohnsonJackson");


    // Constructors

        //- Construct from components
        JohnsonJacksonFrictionalStress(const dictionary& dict);


    //- Destructor
    virtual ~JohnsonJacksonFrictionalStress();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const;

        virtual tmp<volScalarField> muf
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volTensorField& D,
            const dimensionedScalar& phi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::surfaceSlipDisplacementFvPatchField

Description
    fvPatchField corresponding to surfaceSlipDisplacementPointPatchField.
    Is just a slip type since all hard work (projection) is done in the
    pointPatch field.

SourceFiles
    surfaceSlipDisplacementFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSlipDisplacementFvPatchField_H
#define surfaceSlipDisplacementFvPatchField_H

#include "slipFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceSlipDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class surfaceSlipDisplacementFvPatchField
:
    public slipFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("surfaceSlipDisplacement");


    // Constructors

        //- Construct from patch and internal field
        surfaceSlipDisplacementFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        )
        :
            slipFvPatchField<Type>(p, iF)
        {}

        //- Construct from patch, internal field and dictionary
        surfaceSlipDisplacementFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const dictionary& dict
        )
        :
            slipFvPatchField<Type>(p, iF, dict)
        {}

        //- Construct by mapping given surfaceSlipDisplacementFvPatchField
        //  onto a new patch
        surfaceSlipDisplacementFvPatchField
        (
            const surfaceSlipDisplacementFvPatchField<Type>& ptf,
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const fvPatchFieldMapper& mapper
        )
        :
            slipFvPatchField<Type>(ptf, p, iF, mapper)
        {}


        //- Construct as copy
        surfaceSlipDisplacementFvPatchField
        (
            const surfaceSlipDisplacementFvPatchField<Type>& ptf
        )
        :
            slipFvPatchField<Type>(ptf)
        {}

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new surfaceSlipDisplacementFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceSlipDisplacementFvPatchField
        (
            const surfaceSlipDisplacementFvPatchField<Type>& ptf,
            const DimensionedField<Type, volMesh>& iF
        )
        :
            slipFvPatchField<Type>(ptf, iF)
        {}

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new surfaceSlipDisplacementFvPatchField<Type>(*this, iF)
            );
        }

        virtual void write(Ostream& os) const
        {
            slipFvPatchField<Type>::write(os);
            // Only so we can easily postprocess.
            this->writeEntry("value", os);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

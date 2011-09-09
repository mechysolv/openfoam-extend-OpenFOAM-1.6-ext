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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "regionCoupleFvPatchField.H"
#include "magLongDelta.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(iF.name()),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{}


template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(dict.lookup("remoteField")),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{
    if (!isType<regionCoupleFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "regionCoupleFvPatchField<Type>::regionCoupleFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not regionCouple type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const regionCoupleFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_(),
    originalPatchField_(ptf.originalPatchField()),
    curTimeIndex_(ptf.curTimeIndex_)
{
    if (!isType<regionCoupleFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "regionCoupleFvPatchField<Type>::regionCoupleFvPatchField\n"
            "(\n"
            "    const regionCoupleFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const regionCoupleFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    regionCoupleLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(ptf.patch())),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_(),
    originalPatchField_(ptf.originalPatchField()),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a named shadow patch field
template<class Type>
template<class LookupField, class LookupType>
const typename LookupField::PatchFieldType&
regionCoupleFvPatchField<Type>::lookupShadowPatchField
(
    const word& name,
    const LookupField*,
    const LookupType*
) const
{
    // Lookup neighbour field
    const LookupField& shadowField =
        regionCouplePatch_.shadowRegion().
        objectRegistry::lookupObject<LookupField>(name);

    return shadowField.boundaryField()[regionCouplePatch_.shadowIndex()];
}


// Return shadow patch field
template<class Type>
const regionCoupleFvPatchField<Type>&
regionCoupleFvPatchField<Type>::shadowPatchField() const
{
    // Lookup neighbour field
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    return refCast<const regionCoupleFvPatchField<Type> >
    (
        lookupShadowPatchField<GeoField, Type>(remoteFieldName_)
    );
}


// Return neighbour field given internal cell data
template<class Type>
tmp<Field<Type> > regionCoupleFvPatchField<Type>::patchNeighbourField() const
{
    return regionCouplePatch_.interpolate
    (
        shadowPatchField().patchInternalField()
    );
}


// Return neighbour field given internal cell data
template<class Type>
tmp<Field<Type> > regionCoupleFvPatchField<Type>::patchNeighbourField
(
    const word& name
) const
{
    // Lookup neighbour field
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    return regionCouplePatch_.interpolate
    (
        lookupShadowPatchField<GeoField, Type>(name).patchInternalField()
    );
}



template<class Type>
void regionCoupleFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
void regionCoupleFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        originalPatchField_ = *this;
        curTimeIndex_ = this->db().time().timeIndex();
    }


    // Implement weights-based stabilised harmonic interpolation using
    // magnitude of type
    // Algorithm:
    // 1) calculate magnitude of internal field and neighbour field
    // 2) calculate harmonic mean magnitude
    // 3) express harmonic mean magnitude as: mean = w*mOwn + (1 - w)*mNei
    // 4) Based on above, calculate w = (mean - mNei)/(mOwn - mNei)
    // 5) Use weights to interpolate values

    const fvPatch& p = this->patch();
    const scalarField& dc = p.deltaCoeffs();
    const regionCoupleFvPatchField<Type>& spf = this->shadowPatchField();
    const fvPatch& sp = spf.patch();

    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());
    const magLongDelta& smld = magLongDelta::New(sp.boundaryMesh().mesh());

    Field<Type> fOwn =
        this->originalPatchField()
        /(1 - p.weights())/mld.magDelta(p.index());
    Field<Type> fNei = regionCouplePatch_.interpolate
    (
        spf.originalPatchField()
        /(1 - sp.weights())/smld.magDelta(sp.index())
    );

    // Calculate internal weights using field magnitude
    scalarField weights(fOwn.size());

    forAll (weights, faceI)
    {
        scalar mOwn = mag(fOwn[faceI]);
        scalar mNei = mag(fNei[faceI]);

        scalar den = mOwn - mNei;

        if (mag(den) > SMALL)
        {
            scalar mean = mOwn*mNei/(mOwn + mNei)/dc[faceI];
            weights[faceI] = (mean - mNei)/den;
        }
        else
        {
            weights[faceI] = 0.5;
        }
    }

    // Do interpolation
    Field<Type>::operator=(weights*fOwn + (1.0 - weights)*fNei);

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
tmp<Field<Type> > regionCoupleFvPatchField<Type>::snGrad() const
{
    if(regionCouplePatch_.attached())
    {
        return coupledFvPatchField<Type>::snGrad();
    }
    else
    {
        return fvPatchField<Type>::snGrad();
    }
}


// Initialise neighbour processor internal cell data
template<class Type>
void regionCoupleFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes
) const
{
    matrixUpdateBuffer_ = this->patch().patchInternalField(psiInternal);
}


// Return matrix product for coupled boundary
template<class Type>
void regionCoupleFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    scalarField pnf =
        regionCouplePatch_.interpolate
        (
            this->shadowPatchField().matrixUpdateBuffer()
        );

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = regionCouplePatch_.faceCells();

    forAll(fc, elemI)
    {
        result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// Write
template<class Type>
void regionCoupleFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

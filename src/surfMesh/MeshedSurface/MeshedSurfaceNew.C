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

\*---------------------------------------------------------------------------*/

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Face>
Foam::autoPtr< Foam::MeshedSurface<Face> >
Foam::MeshedSurface<Face>::New(const fileName& name, const word& ext)
{
    if (debug)
    {
        Info<< "MeshedSurface::New(const fileName&, const word&) : "
            "constructing MeshedSurface"
            << endl;
    }

    typename fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(ext);

    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        // no direct reader, delegate if possible
        wordHashSet supported = FriendType::readTypes();
        if (supported.found(ext))
        {
            // create indirectly
            autoPtr< MeshedSurface<Face> > surf(new MeshedSurface<Face>);
            surf().transfer(FriendType::New(name, ext)());

            return surf;
        }

        // nothing left to try, issue error
        supported += readTypes();

        FatalErrorIn
        (
            "MeshedSurface<Face>::New(const fileName&, const word&) : "
            "constructing MeshedSurface"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << nl
            << supported
            << exit(FatalError);
    }

    return autoPtr< MeshedSurface<Face> >(cstrIter()(name));
}


template<class Face>
Foam::autoPtr< Foam::MeshedSurface<Face> >
Foam::MeshedSurface<Face>::New(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return New(name, ext);
}

// ************************************************************************* //

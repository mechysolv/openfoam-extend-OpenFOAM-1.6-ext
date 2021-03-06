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

#include "error.H"
#include "sigInt.H"
#include "JobInfo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigInt::oldAction_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigInt::sigIntHandler(int)
{
    // Reset old handling
    if (sigaction(SIGINT, &oldAction_, NULL) < 0)
    {
        FatalErrorIn
        (
            "Foam::sigInt::sigIntHandler()"
        )   << "Cannot reset SIGINT trapping"
            << abort(FatalError);    
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    // Throw signal (to old handler)
    raise(SIGINT);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigInt::sigInt()
{
    oldAction_.sa_handler = NULL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigInt::~sigInt()
{
    // Reset old handling
    if (sigaction(SIGINT, &oldAction_, NULL) < 0)
    {
        FatalErrorIn
        (
            "Foam::sigInt::~sigInt()"
        )   << "Cannot reset SIGINT trapping"
            << abort(FatalError);    
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigInt::set(const bool verbose)
{
    if (oldAction_.sa_handler)
    {
        FatalErrorIn
        (
            "Foam::sigInt::set()"
        )   << "Cannot call sigInt::set() more than once"
            << abort(FatalError);
    }

    struct sigaction newAction;
    newAction.sa_handler = sigIntHandler;
    newAction.sa_flags = SA_NODEFER;
    sigemptyset(&newAction.sa_mask);
    if (sigaction(SIGINT, &newAction, &oldAction_) < 0)
    {
        FatalErrorIn
        (
            "Foam::sigInt::set()"
        )   << "Cannot set SIGINT trapping"
            << abort(FatalError);    
    }
}


// ************************************************************************* //

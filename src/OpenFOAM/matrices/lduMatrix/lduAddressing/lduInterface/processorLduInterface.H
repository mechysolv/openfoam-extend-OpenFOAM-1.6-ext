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
    Foam::processorLduInterface

Description
    An abstract base class for processor coupled interfaces.

SourceFiles
    processorLduInterface.C
    processorLduInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef processorLduInterface_H
#define processorLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorLduInterface Declaration
\*---------------------------------------------------------------------------*/

class processorLduInterface
{
    // Private data

        //- Send buffer.
        //  Only sized and used when compressed or non-blocking comms used.
        mutable List<char> sendBuf_;

        //- Receive buffer.
        //  Only sized and used when compressed or non-blocking comms used.
        mutable List<char> receiveBuf_;

        //- Resize the buffer if required
        void resizeBuf(List<char>& buf, const label size) const;


public:

    //- Runtime type information
    TypeName("processorLduInterface");


    // Constructors

        //- Construct null
        processorLduInterface();


    // Destructor

        virtual ~processorLduInterface();


    // Member Functions

        // Access

            //- Return processor number
            virtual int myProcNo() const = 0;

            //- Return neigbour processor number
            virtual int neighbProcNo() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;


        // Transfer functions

            //- Raw send function
            template<class Type>
            void send
            (
                const Pstream::commsTypes commsType,
                const UList<Type>&
            ) const;

            //- Raw field receive function
            template<class Type>
            void receive
            (
                const Pstream::commsTypes commsType,
                UList<Type>&
            ) const;

            //- Raw field receive function returning field
            template<class Type>
            tmp<Field<Type> > receive
            (
                const Pstream::commsTypes commsType,
                const label size
            ) const;


            //- Raw field send function with data compression
            template<class Type>
            void compressedSend
            (
                const Pstream::commsTypes commsType,
                const UList<Type>&
            ) const;

            //- Raw field receive function with data compression
            template<class Type>
            void compressedReceive
            (
                const Pstream::commsTypes commsType,
                UList<Type>&
            ) const;

            //- Raw field receive function with data compression returning field
            template<class Type>
            tmp<Field<Type> > compressedReceive
            (
                const Pstream::commsTypes commsType,
                const label size
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorLduInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

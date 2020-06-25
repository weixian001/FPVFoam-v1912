/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Contributors/Copyright
    2020 Lim Wei Xian <weixian001@e.ntu.edu.sg> NTUsg

\*---------------------------------------------------------------------------*/

#include "fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
    defineRunTimeSelectionTable(fluidThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo(const fvMesh& mesh, const word& phaseName)
:
    basicThermo(mesh, phaseName),
    Z_                                          //added this section
    (
        IOobject
        (
            "Zmix",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    varZ_
    (
        IOobject
        (
            "varZ",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Chi_
    (
        IOobject
        (
            "chi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Srr_
    (
        IOobject
        (
            "Srr",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -3, -1, 0, 0)
    ),
    Qdt_
    (
        IOobject
        (
            "Qdt",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0, 0)
    )
{}



Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicThermo(mesh, dict, phaseName),
    Z_                                          //added this section
    (
        IOobject
        (
            "Zmix",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    varZ_
    (
        IOobject
        (
            "varZ",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Chi_
    (
        IOobject
        (
            "chi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Srr_
    (
        IOobject
        (
            "Srr",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -3, -1, 0, 0)
    ),
    Qdt_
    (
        IOobject
        (
            "Qdt",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0, 0)
    )
{}


Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    basicThermo(mesh, phaseName, dictionaryName),
    Z_                                          //added this section
    (
        IOobject
        (
            "Zmix",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    varZ_
    (
        IOobject
        (
            "varZ",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Chi_
    (
        IOobject
        (
            "chi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Srr_
    (
        IOobject
        (
            "Srr",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -3, -1, 0, 0)
    ),
    Qdt_
    (
        IOobject
        (
            "Qdt",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -3, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}

//added this section

Foam::volScalarField& Foam::fluidThermo::Z()
{
    return Z_;
}

const Foam::volScalarField& Foam::fluidThermo::Z() const
{
    return Z_;
}

Foam::volScalarField& Foam::fluidThermo::varZ()
{
    return varZ_;
}


const Foam::volScalarField& Foam::fluidThermo::varZ() const
{
    return varZ_;
}

Foam::volScalarField& Foam::fluidThermo::Chi()
{
    return Chi_;
}

const Foam::volScalarField& Foam::fluidThermo::Chi() const
{
    return Chi_;
}

Foam::volScalarField& Foam::fluidThermo::Srr()
{
    return Srr_;
}

const Foam::volScalarField& Foam::fluidThermo::Srr() const
{
    return Srr_;
}

Foam::volScalarField& Foam::fluidThermo::Qdt()
{
    return Qdt_;
}

const Foam::volScalarField& Foam::fluidThermo::Qdt() const
{
    return Qdt_;
}

// ************************************************************************* //

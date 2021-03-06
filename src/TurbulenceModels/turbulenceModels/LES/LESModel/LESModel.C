/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "LESModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::LESModel<BasicTurbulenceModel>::LESModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    LESDict_(this->subOrEmptyDict("LES")),
    turbulence_(LESDict_.get<Switch>("turbulence")),
    printCoeffs_(LESDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(LESDict_.optionalSubDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kMin",
            LESDict_,
            sqr(dimVelocity),
            SMALL
        )
    ),

    epsilonMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilonMin",
            LESDict_,
            kMin_.dimensions()/dimTime,
            SMALL
        )
    ),

    omegaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaMin",
            LESDict_,
            dimless/dimTime,
            SMALL
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            LESDict_
        )
    ),

    CvarZ_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CvarZ",
            //LESDict_,
            coeffDict_,
            1.0
        )
    ),

    Cchi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cchi",
            //LESDict_,
            coeffDict_,
            2.0
        )
    ),

    Sc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sc",
            //LESDict_,
            coeffDict_,
            1.0
        )
    ),

    Sct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sct",
            //LESDict_,
            coeffDict_,
            0.9
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::LESModel<BasicTurbulenceModel>>
Foam::LESModel<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    const IOdictionary modelDict
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    const dictionary& dict = modelDict.subDict("LES");

    const word modelType(dict.get<word>("LESModel"));

    Info<< "Selecting LES turbulence model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "LESModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<LESModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::LESModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        LESDict_ <<= this->subDict("LES");
        LESDict_.readEntry("turbulence", turbulence_);

        coeffDict_ <<= LESDict_.optionalSubDict(type() + "Coeffs");

        delta_().read(LESDict_);

        kMin_.readIfPresent(LESDict_);

	// Turburlent Reacting Flow Coeff
        CvarZ_.readIfPresent(coeffDict_);
        Cchi_.readIfPresent(coeffDict_);
        Sc_.readIfPresent(coeffDict_);
        Sct_.readIfPresent(coeffDict_);

        return true;
    }
	return false;
}


template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::correct()
{
    delta_().correct();
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //

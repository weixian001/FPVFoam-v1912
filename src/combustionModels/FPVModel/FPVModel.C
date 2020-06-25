/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft
    2019 Lim Wei Xian <weixian001@e.ntu.edu.sg> NTUsg

\*---------------------------------------------------------------------------*/

#include "FPVModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"
//#include "physicoChemicalConstants.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
FPVModel<ReactionThermo>::FPVModel
(
//    const word& modelType, const fvMesh& mesh, const word& phaseName
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    //CombThermoType(modelType, mesh, phaseName),
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    solver_(tableSolver(this->mesh(), tables(), Lambda_table())),
    Y_(this->thermo().composition().Y()),
    he_(this->thermo().he()),
    mu_(this->thermo().mu()),
    Srr_(this->thermo().Srr()),			//added             
    Qdt_(this->thermo().Qdt()),
    Z_(this->thermo().Z()),
    //Pv_(this->thermo().Pv()),				//added
    varZ_(this->thermo().varZ()),
    Chi_(this->thermo().Chi()),
    ubIF_(this->mesh().cells().size()),
    ubP_(),
    posIF_(this->mesh().cells().size()),
    posP_(),
    //useScalarDissipation_(this->coeffs().lookup("useScalarDissipation")),
    useProgressVariable_(this->coeffs().lookup("useProgressVariable")),		//added
    useMixtureFractionVariance_(this->coeffs().lookup("useMixtureFractionVariance"))
{
	const polyBoundaryMesh& patches = this->mesh().boundaryMesh();
	int patchSize = 0;
    forAll(patches, patchI)
    {
    	const polyPatch& pp = patches[patchI];
    	if (pp.size() > patchSize) patchSize = pp.size();
    }

    ubP_.setSize(patchSize);
    posP_.setSize(patchSize);
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

//template<class CombThermoType>
template<class ReactionThermo>
FPVModel<ReactionThermo>::~FPVModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class CombThermoType>
template<class ReactionThermo>
hashedWordList FPVModel<ReactionThermo>::tables()
{
	hashedWordList tableNames = this->thermo().composition().species();
	//tableNames.append("chi");
	tableNames.append("HeatRR");
	tableNames.append("mu");
	tableNames.append("he");
	tableNames.append("Srr");		//added
	
	return tableNames;
}

//template<class CombThermoType>
template<class ReactionThermo>
word FPVModel<ReactionThermo>::Lambda_table()
{
        word LambdaNames = "chi_lambda_table";
        return LambdaNames;
}


//template<class CombThermoType>
template<class ReactionThermo>
void FPVModel<ReactionThermo>::correct()
{
    
    // bound the progress variable
    scalar LambdaMax = solver_.maxChi();

    const scalarField& ZCells = Z_.primitiveField();
    const scalarField& varZCells = varZ_.primitiveField();
    const scalarField& chiCells = Chi_.primitiveField();

    scalarField& heCells = he_.primitiveFieldRef();
    scalarField& SrrCells = Srr_.primitiveFieldRef();
    scalarField& QdtCells = Qdt_.primitiveFieldRef();
    scalarField& muCells = mu_.primitiveFieldRef();

    //- Update the species and enthalpy field and reaction rate
    //if(this->active())
    {
       scalarList x(3, 0.0);
       scalarList y(3, 0.0);
       double Zeta;
       double Lambda;

       // Interpolate for internal Field 
       forAll(Y_, i)
       {
    	  scalarField& YCells = Y_[i].primitiveFieldRef();

          forAll(ZCells, cellI)
          {
        	 if (i == 0)
        	 {
        	 Zeta = sqrt(varZCells[cellI]/max(ZCells[cellI]*(1 - ZCells[cellI]), SMALL));
                 y[0] = min(Zeta, 0.99);
                 y[1] = ZCells[cellI];
                 y[2] = max(0.0, chiCells[cellI]);

                 List<int> ubCL_ = solver_.upperBounds_Lambda(y);
                 scalarList posCL_ = solver_.position_Lambda(ubCL_, y);

                 scalarList Lambda_ = solver_.interpolateS(ubCL_, posCL_);
                 Lambda = solver_.upperBoundsLambda_(y[2],Lambda_);

                 //if (useProgressVariable_)   x[0] = min(chiMax, max(0,chiCells[cellI]));
                 if (useProgressVariable_)   x[0] = min(Lambda, LambdaMax);
                 if (useMixtureFractionVariance_) x[1] = min(Zeta, 0.99);
                 x[2] = ZCells[cellI];	

                 ubIF_[cellI] = solver_.upperBounds(x);
                 posIF_[cellI] = solver_.position(ubIF_[cellI], x);
                 
                 QdtCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 4));
            	 muCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 3));
            	 heCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 2));
            	 SrrCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 1));		//added
        	 }

        	 YCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], i);
          }
       }

       // Interpolate for patches
       const volScalarField::Boundary& Chi_Bf = Chi_.boundaryField();
       const volScalarField::Boundary& varZ_Bf = varZ_.boundaryField();
       const volScalarField::Boundary& Z_Bf = Z_.boundaryField();

       volScalarField::Boundary& Qdt_Bf = Qdt_.boundaryFieldRef();
       volScalarField::Boundary& mu_Bf = mu_.boundaryFieldRef();
       volScalarField::Boundary& He_Bf = he_.boundaryFieldRef();
       volScalarField::Boundary& Srr_Bf = Srr_.boundaryFieldRef();

       forAll(he_.boundaryField(), patchi)   
       {
          const fvPatchScalarField& pChi = Chi_Bf[patchi];
          const fvPatchScalarField& pvarZ = varZ_Bf[patchi];
          const fvPatchScalarField& pZ = Z_Bf[patchi];

          fvPatchScalarField& pQdt = Qdt_Bf[patchi];
          fvPatchScalarField& pmu = mu_Bf[patchi];
          fvPatchScalarField& pHe = He_Bf[patchi];
          fvPatchScalarField& pSrr = Srr_Bf[patchi];	//added

          forAll(Y_, i)
          {
		fvPatchScalarField& pY = Y_[i].boundaryFieldRef()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                     Zeta = sqrt(pvarZ[facei]/max(pZ[facei]*(1 - pZ[facei]), SMALL));
                     y[0] = min(Zeta, 0.99);
                     y[1] = pZ[facei];
                     y[2] = max(0.0, pChi[facei]);

                     List<int> ubCL_ = solver_.upperBounds_Lambda(y);
                     scalarList posCL_ = solver_.position_Lambda(ubCL_, y);

                     scalarList Lambda_ = solver_.interpolateS(ubCL_, posCL_);
                     Lambda = solver_.upperBoundsLambda_(y[2],Lambda_);

                     //if (useProgressVariable_)   x[0] = min(chiMax, max(0, pChi[facei]));
                     if (useProgressVariable_)   x[0] = min(Lambda, LambdaMax);
                     if (useMixtureFractionVariance_) x[1] = min(Zeta, 0.99);
                     x[2] = pZ[facei];

                     ubP_[facei] = solver_.upperBounds(x);
                     posP_[facei] = solver_.position(ubP_[facei], x);

                     pQdt[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 4));
                     pmu[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 3));
                     pHe[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 2));
                     pSrr[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 1));			//added
             	 }

            	 pY[facei] = solver_.interpolate(ubP_[facei], posP_[facei], i);
             }
          }
       }

       // Calculate thermodynamic Properties
       this->thermo().correct();
       
       //Info << "Done Calculate thermodynamic Properties" << endl;
    }
}

//template<class CombThermoType>
template<class ReactionThermo>
Switch FPVModel<ReactionThermo>::correctDensity()
{
	return true;
}

//template<class CombThermoType>
template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
FPVModel<ReactionThermo>::R
(
    volScalarField& Y              
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    return tSu;
}

//template<class CombThermoType>
template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FPVModel<ReactionThermo>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tSh;
}

//template<class CombThermoType>
template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FPVModel<ReactionThermo>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tQdot;

/*    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );*/

}

//template<class CombThermoType>
template<class ReactionThermo>
bool FPVModel<ReactionThermo>::read()
{
    if (ThermoCombustion<ReactionThermo>::read())
    {
        this->coeffs().lookup("useProgressVariable") >> useProgressVariable_;
        this->coeffs().lookup("useMixtureFractionVariance") >> useMixtureFractionVariance_;
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

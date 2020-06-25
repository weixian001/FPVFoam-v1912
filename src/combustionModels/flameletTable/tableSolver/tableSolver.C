/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft
    2019 Lim Wei Xian <weixian001@e.ntu.edu.sg> NTUsg

\*---------------------------------------------------------------------------*/


#include "tableSolver.H"

namespace Foam
{
namespace combustionModels
{

defineTypeNameAndDebug(tableSolver, 0);
defineRunTimeSelectionTable(tableSolver, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tableSolver::tableSolver(const fvMesh& mesh, const wordList& tableNames, const word& LambdaNames)
:
   IOdictionary
   (
      IOobject
      (
         "tableProperties",
         mesh.time().constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
      )
   ),
   tableNames_(tableNames),		//tableNames=species+he+Srr+T+rho+mu+psi
   paramNames_(3),
   PreparamNames_(2),
   tables_(tableNames_.size() + 2),
   Lambdatables_(1)
   //chiExt_(this->lookupOrDefault("chiExt", 0.0))
{
    forAll(tableNames_, i)
    {
        tableNames_[i] = tableNames_[i] + "_table";
        tables_.set(i, flameletTable::New(mesh, tableNames_[i])); 
        // return Ptr of flameletTable and store in PtrList<flameletTable> tables_
        //Info << "tables_" << tables_ << endl;
    }

    paramNames_[0] = "Pv_param";
    paramNames_[1] = "Zeta_param";
    paramNames_[2] = "Z_param";

    Lambdatables_.set(0, flameletTable::New(mesh, LambdaNames));
    PreparamNames_[0] = "Zeta_param";
    PreparamNames_[1] = "Z_param";

    forAll(paramNames_, i)
    {
    	params_.append(scalarList(this->lookup(paramNames_[i])));
    }

    forAll(PreparamNames_, i)
    {
        Preparams_.append(scalarList(this->lookup(PreparamNames_[i])));
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tableSolver::~tableSolver()
{}

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

List<int> tableSolver::upperBounds(const scalarList& x) const
{
    // Upper, lower bound for table Interpolation and temp Value for bisection method
	List<int> ub(paramNames_.size(), 0);
    int tlb, tub;
	int newVal = 0;

    // Determine upper bounds and interpolation weights for table interpolation
    // Bisection Method
    for (register label j=0; j<paramNames_.size(); j++)
    {
	   tub = params_[j].size() - 1;
   	   tlb = 0;
   	   while (tub-1 != tlb)
       {
          newVal = (tub + tlb)/2;
          if (x[j] < params_[j][newVal])
     	      tub = newVal;
          else
      	      tlb = newVal;
       }
   	   ub[j] = tub;
    }
   	// Return upper bounds
    return ub;
}

scalarList tableSolver::position(const List<int>& ub, const scalarList& x) const
{
    scalarList pos(paramNames_.size(), 0.0);

    for (register label j=0; j<paramNames_.size(); j++)
    {
        pos[j] = (x[j] - params_[j][ub[j]-1]) /( max( (params_[j][ub[j]] - params_[j][ub[j]-1]), SMALL));
    }

    return pos;
}

scalar tableSolver::interpolate(const List<int>& ub , const scalarList& pos, const label& i) const
{
    //Info << "tables_" << tables_[i] << endl;
    return tables_[i].interpolate(ub, pos);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%% chi-Lamda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List<int> tableSolver::upperBounds_Lambda(const scalarList& y) const
{
        List<int> ub(PreparamNames_.size(), 0);
    int tlb, tub;
        int newVal = 0;

    for (register label j=0; j<PreparamNames_.size(); j++)
    {
           tub = Preparams_[j].size() - 1;
           tlb = 0;
           while (tub-1 != tlb)
       {
          newVal = (tub + tlb)/2;
          if (y[j] < Preparams_[j][newVal])
              tub = newVal;
          else
              tlb = newVal;
       }
           ub[j] = tub;
    }
    return ub;
}

double tableSolver::upperBoundsLambda_(const scalar& y, const scalarList& Lambda_) const
{
        int ub;
    int tlb, tub;
        int newVal = 0;
    if (y < 1.0e-10)
        return params_[0][0];
    else
{
           tub = Lambda_.size() - 1;
           tlb = 0;
           while (tub-1 != tlb)
       {
          newVal = (tub + tlb)/2;
          if (y < Lambda_[newVal])
              tub = newVal;
          else
              tlb = newVal;
       }
           ub = tub;
        double pos = ( y - Lambda_[ub-1]) / ( max( (Lambda_[ub] - Lambda_[ub-1]), SMALL) );
        return params_[0][ub-1]*(1-pos) + params_[0][ub]*pos;
}

}

scalarList tableSolver::position_Lambda(const List<int>& ub, const scalarList& y) const
{
    scalarList pos(PreparamNames_.size(), 0.0);

    for (register label j=0; j<PreparamNames_.size(); j++)
    {
        pos[j] = (y[j] - Preparams_[j][ub[j]-1]) /( max( (Preparams_[j][ub[j]] - Preparams_[j][ub[j]-1]), SMALL) );
    }

    return pos;
}

scalarList tableSolver::interpolateS(List<int>& ub , scalarList& pos) const
{
        return Lambdatables_[0].interpolateS(params_, ub, pos);
}

// %%%%%%%%%%%%%%%%%%%%%%%%% chi-Lambda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scalar tableSolver::maxChi()
{
        /*forAll(paramNames_, i)
        if (paramNames_[i] == "chi_param") return *params_[i].rbegin();
        return 0.0;*/
        return max(params_[0]);
        
}

scalar tableSolver::minChi()
{
        return min(params_[0]);

}

int tableSolver::sizeTableNames() const
{
        return tableNames_.size();
}

} // End combustionModels namespace
} // End Foam namespace

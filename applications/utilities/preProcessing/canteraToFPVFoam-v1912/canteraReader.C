/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    2014 Gabriele Frank <gabriele.frank@unibw.de> Universität der Bundeswehr München
    2019 Lim Wei Xian <weixian001@e.ntu.edu.sg> NTUsg

\*---------------------------------------------------------------------------*/
#include "canteraReader.H"
#include "IFstream.H"
#include "fileName.H"
#include "OFstream.H"
#include "IOList.H"
#include "Gamma.h"
#include <limits>
#include <sstream>

void Foam::canteraReader::read(const fileName& canteraFileName)
{
    //clear the data containers
    singleData_.clear();
    singleData_.setSize(tablesToBeRead_.size());
    tableSorted_.clear();
    coordinates_.clear();
    enthalpyCantera_.clear();

    std::ifstream myinfile(canteraFileName.c_str());    
    yy_buffer_state* bufferPtr(yy_create_buffer(&myinfile, yyBufSize));
    yy_switch_to_buffer(bufferPtr);
    col_iter=0;
    num_columns=0;
    num_lines=0;
    noSecondLine=1;
    while (lex()!=0)
    {}
    yy_delete_buffer(bufferPtr);
}

//void Foam::canteraReader::betaPDFIntegration(const label& numC, const scalar& Zeta)

void Foam::canteraReader::betaPDFIntegration(const label& numPv, const scalar& Zeta)
{
if (Zeta != 0)
{
   List<scalar> Z_(integratedData_[tableNames_["Z"]]);
   List<scalar> varZ_(integratedData_[tableNames_["Z"]].size(), 0.0);
   List<scalar> pdfAlpha(integratedData_[tableNames_["Z"]].size(), 0.0);
   List<scalar> pdfBeta(integratedData_[tableNames_["Z"]].size(), 0.0);
   List<scalar> PDF(integratedData_[tableNames_["Z"]].size(), 0.0);

   for (int i=0;i<integratedData_[tableNames_["Z"]].size();i++)
   {
      varZ_[i] = sqr(Zeta) * (Z_[i]*(1.0 - Z_[i]));
      

      if (varZ_[i] > 1e-7)
      {
          pdfAlpha[i] = Z_[i] * ((Z_[i]*(1.0-Z_[i]))/varZ_[i]-1);
          pdfBeta[i] = (1.0 - Z_[i]) * ((Z_[i]*(1.0-Z_[i]))/varZ_[i]-1);
          
          // Limit alpha and beta but keep their ratio
          if (pdfAlpha[i] > 500)
          {
             pdfBeta[i] = pdfBeta[i]/pdfAlpha[i] * 500;
             pdfAlpha[i] = 500;
          }
          if (pdfBeta[i] > 500)
          {
             pdfAlpha[i] = pdfAlpha[i]/pdfBeta[i] * 500;
             pdfBeta[i] = 500;
          }

          int    gridPoints = 250;
          List<long double> hZ_(gridPoints, 0.0);
          List<long double> helpZ_(gridPoints, 0.0);
          List<long double> delta_(gridPoints, 0.0);

          if ((pdfAlpha[i] > 1) && (pdfBeta[i] > 1))
          {
             // Allocation of Z for PDF integration
             scalar Zmax = 0;
             int    n1 = 0;
             int    n2 = 0;
             PDF.clear();
             PDF.resize(Z_.size(), 0.0);
             
             for (int j=0;j<Z_.size();j++)
             {
                PDF[j] = std::pow(Z_[j],(pdfAlpha[i]-1.0)) * std::pow((1.0 - Z_[j]),(pdfBeta[i]-1.0)) * Gamma(pdfAlpha[i] + pdfBeta[i]) / (min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
                if (PDF[j] > PDF[max(j-1, 0)]) Zmax = Z_[j];
             }
             if(pdfAlpha[i]/pdfBeta[i] <= 0.5)
             {
                n1 = 0.2*gridPoints;
                n2 = 0.8*gridPoints+1;
             }
             else if (pdfAlpha[i]/pdfBeta[i] >= 2)
             {
                n1 = 0.8*gridPoints;
                n2 = 0.2*gridPoints+1;
             }
             else
             {
                n1 = 0.5*gridPoints;
                n2 = 0.5*gridPoints+1;
             }

             //  Allocate Z for 0 < Z < Zmax
             scalar ex1 = 0.9;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));
             

             for (int j=0; j<n1; j++)
             {
                delta_[j] = pow(ex1,j) * delta_[0];
                hZ_[j+1] = hZ_[j] + delta_[j];
             }
             for (int j=1; j<n1; j++)
             {
                hZ_[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 1.1;
             delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
               delta_[j] = pow(ex2,j) * delta_[0];
               helpZ_[j+1] = helpZ_[j] + delta_[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ_[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hZ_[j] /= hZ_[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hZ_.size(), 0.0);

             for (int j=0;j<hZ_.size();j++)
             {
                PDF[j] = (std::pow(hZ_[j],(pdfAlpha[i]-1e0))) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1e0)) * Gamma(pdfAlpha[i] + pdfBeta[i])/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
             }
          }

          else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] > 1))
          {
             // PDF Singularity at Z = 0
             // Allocation of Z for PDF integration
             int    n1 = 0;
             int    n2 = 0;

             if (pdfAlpha[i]/pdfBeta[i] > 0.5)
             {
                scalar Zmax = 0.5;
                scalar ex1 = 1.1;
                n1 = 0.7*gridPoints;
                delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

                // Allocate Z for 0 < Z < Zmax
                for (int j=0; j<n1; j++)
                {
                   delta_[j] = pow(ex1,j) * delta_[0];
                   hZ_[j+1] = hZ_[j] + delta_[j];
                }
                for (int j=1; j<n1; j++)
                {
                   hZ_[j] *= Zmax;
                }

                // Allocate Z for Zmax < Z < 1
                scalar ex2 = 1.1;
                n2 = 0.3*gridPoints+1;
                delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

                for (int j=0; j<n2-1; j++)
                {
                   delta_[j] = pow(ex2,j) * delta_[0];
                   helpZ_[j+1] = helpZ_[j] + delta_[j];
                }
                for (int j=0;j<n2;j++)
                {
                   helpZ_[j] *= (1.0 - Zmax);
                }
                for (int j=0;j<n2-1;j++)
                {
                   hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
                }
             }

             else
	         {
                scalar ex2 = 1.05;
                delta_[0] = (1.0 - ex2)/(1.0 - pow(ex2,(gridPoints-1)));
                for (int j=0; j<gridPoints-1; j++)
                {
                   delta_[j] = pow(ex2,j) * delta_[0];
                   hZ_[j+1] = hZ_[j] + delta_[j];
                }
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hZ_[j] /= hZ_[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hZ_.size(), 0.0);
             for (int j=1;j<hZ_.size();j++)
             {
                PDF[j] = std::pow(hZ_[j],(pdfAlpha[i]-1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
             }
             PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];
          }

          else if ((pdfAlpha[i] > 1) && (pdfBeta[i] <= 1))
          {
          // PDF Singularity at Z = 1
          // Allocation of Z for PDF integration

          int    n1 = 0;
          int    n2 = 0;

          if (pdfAlpha[i]/pdfBeta[i] < 2)
          {
             scalar Zmax = 0.5;
             scalar ex1 = 1.1;
             n1 = 0.3*gridPoints;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             // Allocate Z for 0 < Z < Zmax
             for (int j=0; j<n1; j++)
             {
                delta_[j] = pow(ex1,j) * delta_[0];
                hZ_[j+1] = hZ_[j] + delta_[j];
             }
             for (int j=1; j<n1; j++)
             {
                hZ_[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 0.9;
             n2 = 0.7*gridPoints+1;
             delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
                delta_[j] = pow(ex2,j) * delta_[0];
                helpZ_[j+1] = helpZ_[j] + delta_[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ_[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
             }
          }

          else
          {
             scalar ex1 = 0.95;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(gridPoints-1)));

             for (int j=0; j<gridPoints-1; j++)
             {
                 delta_[j] = pow(ex1,j) * delta_[0];
                 hZ_[j+1] = hZ_[j] + delta_[j];
             }
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hZ_[j] /= hZ_[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hZ_.size(), 0.0);

          for (int j=0;j<hZ_.size()-1;j++)
          {
              PDF[j] = std::pow((hZ_[j]),(pdfAlpha[i]- 1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];
       } 

       else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] <= 1))
       {
          // PDF Singularity at Z = 1 and Z = 0
          // Allocation of Z for PDF integration
          int    n1 = 0;
          int    n2 = 0;
          scalar Zmax = 0.5;
          scalar ex1 = 1.1;
          n1 = 0.5*gridPoints;
          delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

          // Allocate Z for 0 < Z < Zmax
          for (int j=0; j<n1; j++)
          {
              delta_[j] = pow(ex1,j) * delta_[0];
              hZ_[j+1] = hZ_[j] + delta_[j];
          }
          for (int j=1; j<n1; j++)
          {
             hZ_[j] *= Zmax;
          }

          // Allocate Z for Zmax < Z < 1
          scalar ex2 = 0.9;
          n2 = 0.5*gridPoints+1;
          delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

          for (int j=0; j<n2-1; j++)
          {
              delta_[j] = pow(ex2,j) * delta_[0];
              helpZ_[j+1] = helpZ_[j] + delta_[j];
          }
          for (int j=0;j<n2;j++)
          {
              helpZ_[j] *= (1.0 - Zmax);
          }
          for (int j=0;j<n2-1;j++)
          {
             hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hZ_[j] /= hZ_[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hZ_.size(), 0.0);

          for (int j=1;j<hZ_.size()-1;j++)
          {
             PDF[j] = std::pow(hZ_[j],(pdfAlpha[i]- 1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];
          PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];
       }
       
       // Calculate the area of the PDF for scaling
       scalar intPDF = 0;
       for (int j=1;j<hZ_.size();j++)
       {
    	  intPDF += (hZ_[j-1] - hZ_[j]) * (PDF[j-1] + PDF[j])/2;
       }

       // Interpolate singleData entries to the new mixture fraction space
       List<scalar> hY_(gridPoints, 0.0);
       scalar intY = 0;
       
       for (int j=0;j<integratedData_.size();j++)
       {
          hY_ = 0.0;
          hY_[0] = singleData_[j][0];
          hY_[hY_.size()-1] = singleData_[j][singleData_[j].size()-1];
          intY = 0;

          for (int k=1;k<hZ_.size()-1;k++)
          {
             int ubZ = 0;
             for (int l=0;l<Z_.size();l++)
             {
                ubZ = l;
                if (hZ_[k] < Z_[l])
                break;
             }
             int lbZ = ubZ -1;

             // Interpolation to hZ space
             hY_[k] = (singleData_[j][ubZ] - singleData_[j][lbZ])/max(Z_[ubZ] - Z_[lbZ], SMALL) * (hZ_[k] - Z_[lbZ]) + singleData_[j][lbZ];
             // PDF Integration using the trapezoidal rule
             intY += (hZ_[k-1] - hZ_[k]) * (hY_[k-1]*PDF[k-1] + hY_[k]*PDF[k])/(2.0 * intPDF);
          }

          // Special treatment for the boundaries
          intY += (hZ_[hZ_.size()-2] - hZ_[hZ_.size()-1]) * (hY_[hZ_.size()-2]*PDF[hZ_.size()-2] + hY_[hZ_.size()-1]*PDF[hZ_.size()-1])/(2.0 * intPDF);
          if (i != 0 && i != integratedData_[tableNames_["Z"]].size()-1 && tableNames_[j] != "Z")
          {   
          	integratedData_[j][i] = intY; 
          }
       }
     }
   }
 }
}

void Foam::canteraReader::interpolateData()
{
List<scalar> Z_(integratedData_[tableNames_["Z"]]);
List<scalar> hY_(Z_param_.size(), 0.0);

for (int i=0;i<tableNames_.size();i++)
{
   hY_ = 0;
   for (int j=0;j<Z_param_.size();j++)
   {
      int ubZ = 0;
      for (int k=0;k<Z_.size();k++)
      {
         ubZ = k;
         if (Z_[k] > Z_param_[j])
            break;
      }
      int lbZ = ubZ - 1;
      hY_[j] = (integratedData_[i][ubZ] - integratedData_[i][lbZ])/max(Z_[ubZ] - Z_[lbZ], SMALL) * (Z_param_[j] - Z_[lbZ]) + integratedData_[i][lbZ];
    }
    hY_[0] = integratedData_[i][0];
    hY_[Z_param_.size()-1] = integratedData_[i][Z_.size()-1];
    integratedData_[i] = hY_;
}
}

void Foam::canteraReader::calculateEnthalpy()
{
Info << "calculate sensible Enthalpy" << endl;

scalar pstd = 1e5;
label labelT(tableNames_["T"]);
List<scalar> he(singleData_[labelT].size(), 0.0);

for (int i=0;i<singleData_[labelT].size();i++)
{
   for (int j=0; j<thermo.composition().species().size();j++)
   {
      label k = composition.species()[tableNames_[j]];
      he[i] += singleData_[j][i] * composition.Hs(k, pstd, singleData_[labelT][i]);
   }
}
singleData_.append(he);
}

void Foam::canteraReader::calculateZ()
{

if (mixtureFractionDefinition_ == "readFromTable")
{
   Info << "read mixture Fraction from Table" << endl;
}

else
{
   FatalErrorIn("Foam::canteraReader::calculateZ()")
   << "Unknown mixture fraction definition " << nl
   << "Valid mixture fraction definition are :" << nl
   << "readFromTable"  << nl
   << exit(FatalIOError);
}
}

//	added calculation for Progress Variables and the Source Terms

void Foam::canteraReader::calculatePv()
{
   Info << "calculate Progress Variable" << endl;
   
   label labelZ(tableNames_["Z"]);
   List<scalar> Pv(singleData_[labelZ].size(), 0.0);

   for (int i=0;i<singleData_[labelZ].size();i++)
   {
      for (int j=0; j<speciesName_.size();j++)
      {
         label k = tableNames_[speciesName_[j]];
         Pv[i] += singleData_[k][i] * SpeWeightage_[j];
      }
   }
   singleData_.append(Pv); 
}

void Foam::canteraReader::calculateSource()
{
   Info << "calculate Source terms" << endl;
   
   label labelZ(tableNames_["Z"]);
   List<scalar> Srr(singleData_[labelZ].size(), 0.0);

   for (int i=0;i<singleData_[labelZ].size();i++)
   {
      for (int j=0; j<speciesName_.size();j++)
      {
         word speNPortion = "rxn"+speciesName_[j];
         
         label k = tableNames_[speNPortion];
          
         Srr[i] += singleData_[k][i] * SpeWeightage_[j];
      }
   }
   singleData_.append(Srr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::canteraReader::canteraReader(const IOdictionary& canteraDict, const IOdictionary& sPortionDict, rhoReactionThermo& thermo, basicSpecieMixture& composition) //amended
:  composition(composition),
   thermo(thermo),
   Y_(thermo.composition().Y()),
   he_(thermo.he()),         
   p_(thermo.p()),
   tableNames_(thermo.composition().species()),
   tablesToBeRead_(thermo.composition().species()),
   Pv_param_(canteraDict.lookup("Pv_param")),			//amended
   Zeta_param_(canteraDict.lookup("Zeta_param")),
   Z_param_(canteraDict.lookup("Z_param")),
   speciesName_(sPortionDict.lookup("ALL_species")),		//added
   SpeWeightage_(sPortionDict.lookup("Weightages")),		//added
   mixtureFractionDefinition_(canteraDict.lookup("mixtureFractionDefinition")),
   columns_(tableNames_.size()*2+6)
{
   p_ = dimensionedScalar("p", dimPressure, canteraDict.lookup("operatingPressure"));

   //added to read and write the reaction rate of each species
   for (int i=0; i<speciesName_.size(); i++)
   {
	word	rxnSpe = "rxn" + speciesName_[i];	
	tablesToBeRead_.append(rxnSpe);
   	tableNames_.append(rxnSpe);
   }
   
   tablesToBeRead_.append("T");
   tablesToBeRead_.append("Z");
   tablesToBeRead_.append("HeatRR");
   //tablesToBeRead_.append("Cp");
   //tablesToBeRead_.append("psi");
   //tablesToBeRead_.append("rho");
   tablesToBeRead_.append("mu");
   
   tableNames_.append("T");
   tableNames_.append("Z");
   tableNames_.append("HeatRR");
   //tableNames_.append("Cp");
   //tableNames_.append("psi");
   //tableNames_.append("rho");
   tableNames_.append("mu");
   tableNames_.append("he");
   tableNames_.append("chi");
   tableNames_.append("Srr");
  
   fileName canteraFileName(canteraDict.lookup("canteraFileName"));

   //sampledData_ correct size of the Lists
   sampledData_.resize(tableNames_.size());

   for (int i=0; i<sampledData_.size(); i++)
   {
      sampledData_[i].resize(Pv_param_.size());
      for (int j=0; j<sampledData_[i].size();j++)
      {
         sampledData_[i][j].resize(Zeta_param_.size());
         for (int k=0; k<sampledData_[i][j].size();k++)
         {
            sampledData_[i][j][k].resize(Z_param_.size());
         }
      }
   }
   
   newData_.resize(tableNames_.size());

   // newData_ correct size of the Lists
   for (int i=0; i<newData_.size(); i++)
   {
      newData_[i].resize(Zeta_param_.size());
      for (int j=0; j<newData_[i].size();j++)
      {
         newData_[i][j].resize(Z_param_.size());
         for (int k=0; k<newData_[i][j].size();k++)
         {
            newData_[i][j][k].resize(Pv_param_.size());
         }
      }
   }

   // find the correct FileName to read
   for (int numPv=0; numPv<Pv_param_.size();numPv++)	//amended
   {
   	//change chiName to PvName
      std::ostringstream PvName;		//amended
      PvName << Pv_param_[numPv];		//amended

      fileName currentFile = canteraFileName + "_" + PvName.str() + ".csv";	//amended
      Info << "Reading: " << currentFile << endl;
      
      // Read tables
      read(currentFile);
      
      // flip the table if necessary
      List<scalar> hY_(singleData_[tableNames_["Z"]].size(), 0.0);

      if (singleData_[tableNames_["Z"]][singleData_[tableNames_["Z"]].size()-1] <= 0.1)
      {
         for (int i=0; i<singleData_.size(); i++)
         {
            hY_ = 0;
            for (int j=0;j<singleData_[tableNames_["Z"]].size();j++)
            {
               hY_[j] = singleData_[i][singleData_[tableNames_["Z"]].size() - 1 - j];
            } 
            singleData_[i] = hY_;
	 }
      }
      
      // Special treatment for boundaries
      singleData_[tableNames_["Z"]][0] = 0;
      singleData_[tableNames_["Z"]][singleData_[tableNames_["Z"]].size()-1] = 1;
      
      // Calculate hs
      calculateEnthalpy();
      Info << "Done he" << endl;

      // Calculate Z
      calculateZ();
      
      Info << "DOne Z" << endl;
      
      calculatePv();
      
      Info << "Done Pv" << endl;
      
      calculateSource();
      
      Info << "Done Source" << endl;

      for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
      {
         // Write Integrated Data in sampledData
         integratedData_ = singleData_;
         betaPDFIntegration(numPv, Zeta_param_[numZeta]);
         interpolateData();
         for (int k=0; k<integratedData_.size();k++)
         {
            sampledData_[k][numPv][numZeta] = integratedData_[k];
         }
      }
   }

        // added to look up lambda by Z, varZ and Pv
   for (int i=0; i<sampledData_.size(); i++)
   {
      for (int j=0; j<sampledData_[i].size();j++)
      {
         for (int k=0; k<sampledData_[i][j].size();k++)
         {
             for (int m=0; m<sampledData_[i][j][k].size();m++)
             {
                newData_[i][k][m][j] = sampledData_[i][j][k][m];
             }
         }
      }
   }

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::canteraReader::~canteraReader()
{}

// * * * * * * * * * * * * * * * * MemberFunctions  * * * * * * * * * * * * * * * //

Foam::hashedWordList Foam::canteraReader::getNames()
{
	return tableNames_;
}

void	Foam::canteraReader::write(const int& i,Foam::IOdictionary& dictionary, Foam::OFstream& output)
{
	word dictionaryName=tableNames_[i]+"_table";
	List<List<List<scalar> > >lists=sampledData_[i];
	dictionary.set(dictionaryName,lists);
	dictionary.writeHeader(output);
	output<<dictionaryName<<lists<<";";
}

void    Foam::canteraReader::writechi(Foam::IOdictionary& dictionaryChi, Foam::OFstream& output)
{
        word ProgressVariableName=tableNames_[tableNames_.size()-2]+"_lambda_table";  //added
        List<List<List<scalar> > >chiArray=newData_[tableNames_.size()-2];         //added
        dictionaryChi.set(ProgressVariableName,chiArray);
        dictionaryChi.writeHeader(output);
        output<<ProgressVariableName<<chiArray<<";";
}


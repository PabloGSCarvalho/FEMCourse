/***************************************************************************
 *   Copyright (C) 2005 by Philippe R. B. Devloo                           *
 *   phil@fec.unicamp.br                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "tanalysis.h"
#include "TMatrix.h"
#include <math.h>
#include "TVec.h"
#include "tpanic.h"

TAnalysis::TAnalysis(TMalha *malha) : fMalha(malha), fSolution()
{
}


TAnalysis::~TAnalysis()
{
}

// Builds the global stiffness matrix and right hand side

void TAnalysis::Run()
{
    TMalha *malha=this->fMalha;
    TVec<TNo> nos=malha->getNodeVec();
    double nnodes=nos.Size();
    TMatrix K(nnodes,nnodes);
    K.Zero();
    TMatrix F(nnodes,1);
    F.Zero();
    
    Assemble(K, F);
    
    K.Print();
    F.Print();
    
    K.Solve_LU(F);
    fSolution=F;
    fSolution.Print();
}


// assembly method

void TAnalysis::Assemble(TMatrix &stiff, TMatrix &rhs)
{
    //Inicializar as matrizes
    
    TMalha *malha=this->fMalha;
    double nel=malha->getElementVec().Size();
    double nnodes=malha->getNodeVec().Size();
    
    stiff.Resize(nnodes, nnodes);
    rhs.Resize(nnodes, 1);
    stiff.Zero();
    rhs.Zero();
    
    for (int el=0; el<nel; el++) {
        
        
        TMalha *malha=this->fMalha;
        TElemento *cel=malha->getElement(el);
        TVec<int> nos=cel->getNodeVec();
        
        TMatrix EK(nos.Size(),nos.Size()),EF(nos.Size(),1);
        
        EF.Zero();
        EK.Zero();
        
        cel->CalcStiff(*malha, EK, EF);
        
        //EK.Print();
        
        for (int i=0; i<nos.Size(); i++) {
            rhs(nos[i],0)+=EF(i,0);
            for (int j=0; j<nos.Size(); j++) {
                stiff(nos[i],nos[j])+=EK(i,j);
            }
        }
        
        //stiff.Print();
    }
    
    
}



// computes the error in energy and l2 norm


void TAnalysis::Error(void (*exact) (TVec<double> &x, double &val, TVec<double> &deriv), double &energy, double &l2)
{
    
    TMalha *malha=this->fMalha;
    double nel=malha->getElementVec().Size();
    
    double erroh1=0.,errol2=0.;
    
    for (int el=0; el<nel; el++) {
       
        TElemento *cel=malha->getElement(el);
        if (cel->GetMatid() <= 0) {
            continue;
        }
            
        cel->Error(fSolution, *fMalha,exact, erroh1, errol2);
        
        energy+=erroh1;
        l2+=errol2;
    }
    
    energy = sqrt(energy);
    l2 = sqrt(l2);    
}

// computes the global solution
void TAnalysis::uh(std::string &file_name){
    
    
    std::ofstream myfile;
    myfile.open(file_name.c_str());
    
    TMalha *malha=this->fMalha;
    double nel=malha->getElementVec().Size();
    int np = 3;
    TMatrix uh;
    TMatrix duhdx;
    for (int el=0; el<nel; el++) {
        
        TElemento *cel=malha->getElement(el);
        if (cel->GetMatid() <= 0) {
            continue;
        }
        
        cel->uhe(fSolution, *malha, uh, duhdx);
        for (int i = 0 ; i< np; i++) {
            myfile << uh(i,0) << " " << uh(i,1) << " " << uh(i,0) << " " << duhdx(i,1) << std::endl;
        }
    }
    
    
    myfile.close();
    
}
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
#ifndef TANALYSIS_H
#define TANALYSIS_H

#include "tmalha.h"
#include "telemento.h"
#include "TMatrix.h"
/*
This class implements a linear analysis for a finite element simulation

@author Philippe R. B. Devloo
*/
class TAnalysis{
public:
    TAnalysis(TMalha *malha);

    ~TAnalysis();

    // Builds the global stiffness matrix and right hand side
    void Run();
    
    // computes the error in energy and l2 norm
    void Error(void (*exact) (TVec<double> &x, double &val, TVec<double> &deriv), double &energy, double &l2);
    
    // computes the global solution
    void uh(std::string &file_name);
    
    void Assemble(TMatrix &stiff, TMatrix &rhs);
 
private:

    // The mesh for which an analysis is performed
    TMalha *fMalha;
    
    // The solution of the problem
    TMatrix fSolution;

    // Pointer to the function which computes the exact solution
//    void (*fFunction) (TVec<double> &x, double &val, TVec<double> &deriv);
    
        
};

#endif

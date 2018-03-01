/***************************************************************************
 *   Copyright (C) 2005 by Philippe R. B. Devloo                           *
 *   phil@corona                                                           *
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
#include "tmaterialbc.h"
#include "TMatrix.h"
#include "TVec.h"
#include "tpanic.h"

static double BigNumber = 1.e12;

TMaterialBC::    TMaterialBC(int id, int bctype, double contrstiff, double contrrhs) : TMaterial(id), fContrStiff(contrstiff), fContrRhs(contrrhs), fBCType(bctype)
{
}


TMaterialBC::~TMaterialBC()
{
}


void TMaterialBC::Print(std::ostream& out) const
{
    TMaterial::Print(out);
    out << "BC type " << fBCType;
    out << " BC values stiff " << fContrStiff << " rhs " << fContrRhs << std::endl;
}

void TMaterialBC::Contribute (double  weight,
                              TVec<double> & philVal,
                              TMatrix & dphi,TMatrix & elementK,
                              TMatrix & elementF) const
{
    
    
    if (fBCType==0) {
    
        //Dirichlet
        for (int i=0; i<philVal.Size(); i++) {
            for (int j=0; j<philVal.Size(); j++) {
                elementK(i,j)+=philVal[i]*philVal[j]*weight*BigNumber;
            }
            elementF(i,0)+=philVal[i]*weight*BigNumber*fContrRhs;
            
        }
        
        
    }else if (fBCType==1){
        
        //Neumman
        for (int i=0; i<philVal.Size(); i++) {
            
            elementF(i,0)+=philVal[i]*weight*fContrRhs;
        }
        
    
    }else if (fBCType==2){
        
        //Mista
        for (int i=0; i<philVal.Size(); i++) {
            for (int j=0; j<philVal.Size(); j++) {
                elementK(i,j)+=philVal[i]*philVal[j]*weight*fContrStiff;
            }
            elementF(i,0)+=philVal[i]*weight*fContrRhs;
            
        }
        
    
    }
    
    
    
    
}


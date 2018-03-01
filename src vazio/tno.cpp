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
#include "tno.h"
#include "tpanic.h"


TNo::TNo()
{
    fCoordenadas[0]=0;
    fCoordenadas[1]=0;
}


TNo::~TNo()
{
}

TNo::TNo (TVec<double> &coordenadas)
{
    if (coordenadas.Size()!=2) {
        DebugStop();
    }
    
    fCoordenadas[0]=coordenadas[0];
    fCoordenadas[1]=coordenadas[1];
    
}
/*
 TVec<double>::iterator it;
 for(it=coordenadas.begin(); it!=coordenadas.end() && i<2; it++,i++)
 {
 fCoordenadas[i] = (*it);
 }
 */

TNo::TNo (const TNo &copy)
{
    fCoordenadas[0]=copy.fCoordenadas[0];
    fCoordenadas[1]=copy.fCoordenadas[1];
    
}


void TNo::setData(TVec<double> &coordenadas)
{
    if (coordenadas.Size()!=2) {
        DebugStop();
    }
    
    fCoordenadas[0]=coordenadas[0];
    fCoordenadas[1]=coordenadas[1];
    
}

double &TNo::Co(int i)
{
    if (i!=0&&i!=1) {
        DebugStop();
    }
    
    return fCoordenadas[i];
}



void TNo::Print(std::ostream &out) const
{
    out << "No: coordenadas ( ";
    int i;
    out << fCoordenadas[0];
    for (i=1;i<2;i++)
    {
        out << " , " << fCoordenadas[i];
    }
    
    out  << " );" << std::endl;
}

bool TNo::operator==(const TNo &other) const
{
    DebugStop();
    return false;
}


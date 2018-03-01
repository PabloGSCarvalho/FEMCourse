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
#include "tmalha.h"
#include "telemento.h"
#include "tmaterial.h"
#include "telemento1d.h"
#include "telemento0d.h"
#include "tmaterialbc.h"
#include "tmaterial1d.h"
#include "TMatrix.h"
#include "tpanic.h"

TMalha::TMalha() : fNodes(), fElements(), fMaterials()
{
}


TMalha::~TMalha()
{
    
}


int TMalha::insertNode(TNo &node)
{
    int result = fNodes.Size();
    fNodes.Resize(result+1);
    fNodes[result]=node;
    return result;
}

int TMalha::insertMaterial(TMaterial *mat)
{
    int id=mat->Id();
    fMaterials[id]=mat;
    
    return id;
}

int TMalha::insertElement (TElemento *el)
{
    int result = fElements.Size();
    fElements.Resize(result+1);
    fElements[result]=el;

    return result;
}

TVec <TNo> & TMalha::getNodeVec()
{
    return fNodes;
}

std::map <int, TMaterial *> & TMalha::getMaterialMap()
{
    return fMaterials;
}

TVec <TElemento *> & TMalha::getElementVec()
{
    return fElements;
}


TNo &TMalha::getNode(int id)
{
    return fNodes[id];
}


TMaterial * TMalha::getMaterial (int i)
{
    return fMaterials[i];
}


TElemento * TMalha::getElement (int i)
{
    return fElements[i];
}


void TMalha::Print(std::ostream &out)
{
    int i;
    out << "Impress‹o da malha" << std::endl;
    out << "Vetor de Nos - número de nos = " << fNodes.Size() << std::endl;
    for (i=0;i<fNodes.Size();i++)
    {
        out << "Ind " << i << ' ';
        fNodes[i].Print(out);
    }
    out << "Vetor de Elementos - numero de elementos = " << fElements.Size() << std::endl;
    for (i=0; i<this->fElements.Size(); i++)
    {
        out << "Ind " << i << ' ';
        fElements[i]->Print(out);
    }
    out << "Vetor de Materiais - numero de materiais = " << fMaterials.size() << std::endl;
    std::map<int, TMaterial *>::iterator it;
    for (it=fMaterials.begin(); it != this->fMaterials.end(); it++)
    {
        (*it).second->Print(out);
    }
}


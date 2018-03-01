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
#include "telemento.h"
#include "tno.h"
#include "tmalha.h"
#include "TMatrix.h"
#include "tpanic.h"

TElemento::TElemento() : fNodes(), fMaterialId(-1), fPorder(-1)
{
    TMatrix a(4,4,0.);
    
}

/**
 *  Elemento criado com
 * @matid : indice do material
 * @order : ordem de interpolacao
 * @nodes : indices dos nos
 * nota que existe correspondencia entre o numero de nos e a ordem do elemento
 */

TElemento::TElemento(int matid, int order, TVec<int> &nodes) :
fNodes(nodes), fMaterialId(matid), fPorder(order)
{
}


TElemento::~TElemento()
{
}

/**
 * Cálcula o valor da função de forma em um dado ponto
 * @param pt [in] ponto onde se quer calcular o valor das funções de forma
 * @param phiValue [out] valor de cada função de forma do elemento no ponto
 */

void TElemento::Shape1d (int order, TVec<double> & pt,TVec<double> & phi, TMatrix& dphi)
{
    
    
    phi.Resize(order+1);
    dphi.Resize(1, order+1);
    for (int i=0; i<order+1; i++) {
        phi[i]=1;
        dphi(0,i)=0;
    }
    
    
    
    for (int i=0; i<order+1; i++) {
        double epsi=-1.+i*2./order;
        double axdphi;
        
        for (int j=0; j<order+1; j++) {
            
            
            if (i!=j) {
                double epsj=-1.+j*2./order;
                
                phi[i]*=(pt[0]-epsj)/(epsi-epsj);
                
                axdphi=1/(epsi-epsj);
                
                for (int k=0; k<order+1; k++) {
                    if (k!=i&&k!=j) {
                        epsj=-1.+k*2./order;
                        axdphi*=(pt[0]-epsj)/(epsi-epsj);
                    }
                    
                }
                
                dphi(0,i)+=axdphi;

                
            }
            
            
        }
        
        
        
    }

}

int TElemento::main()
{
    TVec<double> pt(1), phi(3);
    TMatrix dphi(1,3);
    pt[0] = 0.5;
    Shape1d(3,pt,phi,dphi);
    return 0;
}



/*!
 \fn TElemento::Print(std::ostream &out = cout)
 */

void TElemento::Print(std::ostream &out)
{
    out << "ElType " << TypeName(getElType()) << " matid " << fMaterialId << " order " << fPorder << " nodes ";
    int i;
    for(i=0; i<fNodes.Size(); i++) out << fNodes[i] << ' ';
    out << std::endl;
}



/*
 \fn TElemento::TypeName(MElementType type)
 */

std::string TElemento::TypeName(MElementType type)
{
    std::string result;
    switch (type)
    {
        case EPoint:
            result = "Point";
            break;
        case ELinear:
            result = "Linear";
            break;
        case EQuadrilateral:
            result = "Quadrilateral";
            break;
        case  ETriangle:
            result = "Triangle";
            break;
        default:
            result = "Unknown";
    }
    return result;
}



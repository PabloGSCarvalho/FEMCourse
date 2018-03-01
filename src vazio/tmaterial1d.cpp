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
#include "tmaterial1d.h"
#include "TMatrix.h"
#include "tpanic.h"

TMaterial1d::TMaterial1d(int id, double K, double C, double B, double F) : TMaterial(id), fK(K), fC(C), fB(B), fF(F)
{
}



TMaterial1d::~TMaterial1d()
{
}


void TMaterial1d::Print(std::ostream& out) const
{
    TMaterial::Print(out);
    out << "Coeficient values K " << fK << " C " << fC << " B " << fB << " F " << fF << std::endl;
}

/**
 * Calcula o valor da contribuição da equação variacional no ponto dado
 * na matriz de rigidez do elemento e no vetor de carga
 * @param pt [in]: ponto de integração de Gauss
 * @param weight [in]: peso de integração
 * @param phiVal [in] : valor da função teste no ponto dado
 * @param dphi [in] : valor das derivadas da função de forma no ponto de integração
 * @param elementK [inout]: matriz de rigidez do elemento
 * @param elementF [inout]: vetor de carga do elemento
 */
void TMaterial1d::Contribute (double  weight,
                              TVec<double> & philVal,
                              TMatrix & dphi,TMatrix & elementK,
                              TMatrix & elementF) const
{
    
    for (int i=0; i<philVal.Size(); i++) {
        for (int j=0; j<philVal.Size(); j++) {
            elementK(i,j)+=fK*dphi(0,i)*dphi(0,j)*weight+fB*philVal[i]*philVal[j]*weight+fC*philVal[i]*dphi(0,j)*weight;
        }
        elementF(i,0)+=weight*philVal[i]*fF;
    }
    
    //elementK.Print("Pos Contribute");
}

/**
 * Calcula a contribuicao para o erro da solucao
 * @param weight [in] peso do ponto de integracao
 * @param sol [in] valor da solucao
 * @param deriv [in] valor da derivada da solucao
 * @param function [in] ponteiro para funcao que calcula o valor exato
 * @param energy [in/out] contribuicao para norma da energia
 * @param l2 [in/out] contribuicao para norma em L2
 *
 */
void TMaterial1d::ContributeErrorSquare(TVec<double> &x, double weight, double sol, TVec<double> &deriv,
                                        void (function)(TVec<double>& x, double &val, TVec<double>&der), double &energy, double &l2)
{
    double val;
    TVec<double> deriuexact(1);
    function(x,val,deriuexact);
    l2+=weight*(val-sol)*(val-sol);
    energy+=weight*(deriuexact[0]-deriv[0])*(deriuexact[0]-deriv[0]);
}

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
#ifndef TELEMENTO0D_H
#define TELEMENTO0D_H

#include <telemento.h>
#include "TVec.h"

/**
Esta classe implementa um elemento ponto

@author Philippe R. B. Devloo
*/
class TElemento0d : public TElemento
{
public:
    TElemento0d();
    
    TElemento0d(int matid, int order, TVec< int >& nodes);

    ~TElemento0d();

    static std::string TypeName(MElementType type);
    virtual void CalcStiff(TMalha &malha, TMatrix &stiff, TMatrix &rhs);
    virtual void Jacobian(TVec< double > &point, TMatrix& jacobian, TMatrix& jacinv, double &detjac, TMalha& malha);
    
    // retorno o tipo de elemento
    MElementType getElType();
    
    /**
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
    virtual void Shape(TVec<double> &point, TVec<double> &phi, TMatrix &dphi);
    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
virtual void Error(TMatrix &solution, TMalha &malha, void (*f)(TVec<double> &,double &, TVec<double> &), double &energy, double &l2)
     {
     	energy = 0.;
	l2 = 0.;
     }
    
    /**
     * Calcula o a solucao e o gradiente do elemento
     * @param solution vetor de coeficientes multiplicadores (dof)
     * @param malha espaco de aproximacao
     * @param uhe combinacao linear de alpha_{i} phi_{i}
     * @param duhedx combinacao linear de alpha_{i} dphi_{i}
     */
    virtual void uhe(TMatrix &solution, TMalha &malha, TMatrix &uhe, TMatrix &duhedx);


};

#endif

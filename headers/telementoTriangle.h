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
#ifndef TELEMENTOTRIANGLE_H
#define TELEMENTOTRIANGLE_H

#include <telemento.h>
#include "TVec.h"

/**
@author Philippe R. B. Devloo
*/
class TElementoTriangle : public TElemento
{
public:
    TElementoTriangle(int matid, int order, TVec< int >& nodes);
    
    TElementoTriangle();

    ~TElementoTriangle();

    virtual MElementType getElType();
    virtual void CalcStiff(TMalha &malha, TMatrix& stiff, TMatrix& rhs);
    virtual void Jacobian(TVec<double> &point, TMatrix& jacobian, TMatrix& jacinv, double &detjac, TMalha& malha);
    
    
    /*
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
virtual void Error(TMatrix &solution, TMalha &malha, void (f)(TVec<double> &,double &, TVec<double> &), double &energy, double &l2);

    virtual void uhe(TMatrix &solution, TMalha &malha, TMatrix &uhe, TMatrix &duhedx);

};

#endif

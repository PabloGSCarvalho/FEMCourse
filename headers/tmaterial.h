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
#ifndef TMATERIAL_H
#define TMATERIAL_H

#include <iostream>
#include "TVec.h"

class TMatrix;

/**
Implementa a interface típica de um material


@author Philippe R. B. Devloo
*/
class TMaterial{
public:
    /**
     * Construtor vazio
     */
    TMaterial(int id);
    
    /**
     * Construtor de copia
     */
    TMaterial (const TMaterial & copy);

    /**
     * Destrutor
     */
    virtual ~TMaterial();
    
public: //funções de cálculo

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
  virtual void Contribute (double  weight,
                           TVec<double> & philVal,
                           TMatrix & dphi,TMatrix & elementK,
                           TMatrix & elementF) const = 0; 

   /**
    * Calcula a contribuicao para o erro da solucao
    * @param x [in] localizacao do ponto
    * @param weight [in] peso do ponto de integracao
    * @param sol [in] valor da solucao
    * @param deriv [in] valor da derivada da solucao
    * @param function [in] ponteiro para funcao que calcula o valor exato
    * @param energy [in/out] contribuicao para norma da energia
    * @param l2 [in/out] contribuicao para norma em L2
    *
    */
   virtual void ContributeErrorSquare(TVec<double> &x, double weight, double sol, TVec<double> &deriv,
   	void (function)(TVec<double>& x, double &val, TVec<double>&der), double &energy, double &l2) = 0;
	
public://Diversos
  /**
   * Imprime os dados da classe
   */
  virtual void Print(std::ostream &out = std::cout) const;
  
    static double GetBigNumber()
    {
        return gBigNumber;
    }
    
  /**
   * retorna o identificador do material
   */
   int Id() 
   {
    return fId;
   }
  
protected:
  /**
   * Identificador do material
   */
   int fId;
    
   ///valor do BigNumber a ser adotado
    static const double gBigNumber;
   
};

#endif

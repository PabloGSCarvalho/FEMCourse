//
//  TMaterial2D.hpp
//  FemSC
//
//  Created by Philippe Devloo on 11/12/15.
//
//

#ifndef TMATERIAL2D_H
#define TMATERIAL2D_H

#include <stdio.h>

#include "TVec.h"

#include <tmaterial.h>

/**
This class implements the variational formulation of a 1D partial differential equation

@author Philippe R. B. Devloo
*/
class TMaterial2d : public TMaterial
{
public:

	TMaterial2d(int id, double K, double C[2], double B, double F);

	~TMaterial2d();

	/**
	* Calcula o valor da contribui햫o da equa햫o variacional no ponto dado
	* na matriz de rigidez do elemento e no vetor de carga
	* @param pt [in]: ponto de integra햫o de Gauss
	* @param weight [in]: peso de integra햫o
	* @param phiVal [in] : valor da fun햫o teste no ponto dado
	* @param dphix [in] : valor das derivadas da fun햫o de forma no ponto de integra햫o
	* @param elementK [inout]: matriz de rigidez do elemento
	* @param elementF [inout]: vetor de carga do elemento
	*/
	virtual void Contribute(double  weight,
		TVec<double> & philVal,
		TMatrix & dphix, TMatrix & elementK,
		TMatrix & elementF) const;
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
		void (function)(TVec<double>& x, double &val, TVec<double>&der), double &energy, double &l2);

	virtual void Print(std::ostream& out) const;

protected:

	/**
	* Definition of the differential equation coeficients according to the book of Becker, Carey and Oden
	*/
	double fK, fC[2], fB, fF;

};


#endif /* TMaterial2d_hpp */


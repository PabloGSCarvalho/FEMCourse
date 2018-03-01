//
//  TMaterial2D.cpp
//  FemSC
//
//  Created by Philippe Devloo on 11/12/15.
//
//

#include "TMaterial2d.h"

#include "TMatrix.h"

TMaterial2d::TMaterial2d(int id, double K, double C[2], double B, double F) : TMaterial(id), fK(K), fB(B), fF(F)
{
    fC[0] = C[0];
    fC[1] = C[1];
}



TMaterial2d::~TMaterial2d()
{
}


void TMaterial2d::Print(std::ostream& out) const
{
    TMaterial::Print(out);
    out << "Coeficient values K " << fK << " C " << fC[0] << " " << fC[1] << " B " << fB << " F " << fF << std::endl;
}


//F -> LoadVec coef;
//B -> ''

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

void TMaterial2d::Contribute (double  weight,
                              TVec<double> & philVal,
                              TMatrix & dphix,TMatrix & elementK /*stiff*/,
                              TMatrix & elementF /*rhs*/) const
{
    int i, j, nshape;
    nshape = philVal.Size();
    //  std::cout<<"valor das funcoes de forma no ponto de integracao:"<<std::endl;
    //  philVal.Print();
    //  elementK.Print("K local antes do contribute");
    //  elementF.Print("F local antes do contribute");
    for(i=0; i<nshape; i++)
    {
        elementF(i,0) += weight*philVal[i]*fF;
        for(j=0; j<nshape; j++)
        {
            elementK(i,j) += (dphix(0,i)*dphix(0,j)+dphix(1,i)*dphix(1,j))*fK*weight+
            weight*philVal[i]*(dphix(0,j)*fC[0]+dphix(1,j)*fC[1])+
            philVal[i]*philVal[j]*fB*weight;
        }
    }
    
    //elementK.Print("K local depois do contribute");
    //  elementF.Print("F local depois do contribute");
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
void TMaterial2d::ContributeErrorSquare(TVec<double> &x, double weight, double sol, TVec<double> &deriv,
                                        void (function)(TVec<double>& x, double &val, TVec<double>&der), double &energy, double &l2)
{
    double solexact=0.;
    int dimension = deriv.Size();
    TVec<double> derivsolexact(dimension);
    for (int i=0; i<dimension; i++) {
        derivsolexact[i]=0.0;
    }
    
    function(x,solexact,derivsolexact);
    int id;
    
    for(id=0; id<dimension; id++)
    {
        energy += weight*(deriv[id]-derivsolexact[id])*(deriv[id]-derivsolexact[id]);
    }
    energy += weight*(sol-solexact)*(sol-solexact);
    l2 += weight*(sol-solexact)*(sol-solexact);
    
}

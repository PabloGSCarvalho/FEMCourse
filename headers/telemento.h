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
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *t Wt WIIt WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GN
    
    TNos * getNode(i);
U General Public License for more details.          t WI                *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TELEMENTO_H
#define TELEMENTO_H

class TMalha;

//class TPZMatrix;
//class TMatrix;

#include "TVec.h"
#include "TMatrix.h"
#include <iostream>
#include <string>
/**
Classe virtual que implementa a interface básica de um elemento

@author Philippe R. B. Devloo
*/
/**
 * Definição dos tipos de elementos
 */
 
enum MElementType { EPoint, ELinear, ETriangle, EQuadrilateral };

class TElemento{
public:
    /**
     * Construtor vazio - inicializar com valores que indiquem falta de 
     * inicialização
     */
    TElemento();
    
    /**
     *  Elemento criado com
     * @matid : indice do material
     * @order : ordem de interpolacao
     * @nodes : indices dos nos
     * nota que existe correspondencia entre o numero de nos e a ordem do elemento
     */
     TElemento(int matid, int order, TVec<int> &nodes);
    
    /**
     * Destrutor padrão
     */
    ~TElemento();
    
    /**
     * Metodo para acesso aos nos do elemento
     */
     const TVec<int> &getNodeVec() const
     {
       return fNodes;
     }
    
public: //Acesso a leitura dos dados da classe
    /**
     * Retorna o tipo de elemento
     */
    virtual MElementType getElType()=0;

public: //Métodos de cálculo
    
    static int main();
    
public:

    /**
     * Cálcula o valor da função de forma em um dado ponto
     * @param pt [in] ponto onde se quer calcular o valor das funções de forma
     * @param phiValue [out] valor de cada função de forma do elemento no ponto
     */
    static void Shape1d (int order, TVec<double> & pt,TVec<double> & phi, TMatrix &dphi);
    
    /**
     * Calcula a matriz de rigidez local
     */
    virtual void CalcStiff (TMalha &malha, TMatrix &stiff, TMatrix &rhs) = 0;
    
    /**
     * Calculo do jacobiano
     * @point : ponto em qual calculamos o jacobiano do mapeamento
     * @jacobian : matriz jacobiana (a dimensao depende da dimensao do elemento
     * @jacinv : inverso da matriz jacobiana
     * @malha : objeto malha necessaria para relacionar os indices dos nos com os nos reais
     * 
     */
     virtual void Jacobian(TVec<double> &point, TMatrix &jacobian, TMatrix &jacinv, double &detjac, TMalha &malha) = 0;
     
     /**
      * return the string corresponding to the type
      */
    static std::string TypeName(MElementType type);
    
    /**
     * Print the element
     */
    void Print(std::ostream &out = std::cout);
    
    /**
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
    virtual void Shape(TVec<double> &point, TVec<double> &phi, TMatrix &dphi) = 0;
    
    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
virtual void Error(TMatrix &solution, TMalha &malha, void (*f)(TVec<double> &,double &, TVec<double> &), double &energy, double &l2) = 0;
    
    /**
     * Calcula o a solucao e o gradiente do elemento
     * @param solution vetor de coeficientes multiplicadores (dof)
     * @param malha espaco de aproximacao
     * @param uhe combinacao linear de alpha_{i} phi_{i}
     * @param duhedx combinacao linear de alpha_{i} dphi_{i}
     */
    virtual void uhe(TMatrix &solution, TMalha &malha, TMatrix &uhe, TMatrix &duhedx) = 0 ;
    
    /**
     * Retorna o indice material associado ao elemento
     */
    
    int GetMatid(){
        return fMaterialId;
    }
    
protected:
    
    /**
     * Vetor de ponteiros para os nós. Os nós são armazenados
     * na malha
     */
    TVec<int> fNodes;
    
    /**
     * Referência para o índice do material na malha
     */
    int fMaterialId;
    
    /**
     * Ordem de interpolação p
     */
    int fPorder;
    
};

#endif

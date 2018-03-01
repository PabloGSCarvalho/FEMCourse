//  TVecNum.cpp
//
//  Created by Eduardo on 03/03/15.
//  Copyright (c) 2015 Eduardo. All rights reserved.
//

#include <iostream>
#include "TVecNum.h"
#include <stdio.h>
#include <math.h>
#include "tpanic.h"

using std::cout;
using std::endl;
using std::cin;

//construtor padrao de TVec
template<class T>
TVecNum<T>::TVecNum(): TVec<T>()
{
}

//fvector recebe o tamanho (fsize) em alocacao dinamica
//todos os elementos do vetor recebem zero
template<class T>
TVecNum<T>::TVecNum(int size) : TVec<T>(size)
{
    for (int i=0; i<size; i++) {
        this->fvector[i] = 0;
    }
}

//fvector recebe o tamanho (fsize) em alocacao dinamica
//todos os elementos do vetor recebem zero
template<class T>
TVecNum<T>::TVecNum(int size, T val) : TVec<T>(size)
{
    for (int i=0; i<size; i++) {
        this->fvector[i] = val;
    }
}

///construtor copia
template<class T>
TVecNum<T>::TVecNum(const TVecNum &cp) : TVec<T>(cp)
{
}

///Destrutor dos vetores criados
template<class T>
TVecNum<T>::~TVecNum(void)
{
}


///cria um ponteiro "*copia" que recebe os valores do fvector, modifica o tamanho do vetor e "*copia" devolve os valores para o vetor de novo tamanho;
template<class T>
void TVecNum<T>::Resize (int nsize)
{
    int prevsize = this->fsize;
    TVec<T>::Resize(nsize);
    for (int i=prevsize; i<nsize; i++) {
        this->fvector[i] = 0;
    }
    ///std::cout << "Expansao/Trucamento de veto:";
}

///funcao de ordenacao de vetores, trabalha em conjunto com a funcao Swap
template<class T>
void TVecNum<T>::OrdVec ()
{
    if (this->fvector == NULL) {
        DebugStop();
    }
    
    ///std::cout << "Ordenacao de vetor:";
    for (int i=0; i<this->fsize; i++) {
        
        int menor = i; //menor indice
        
        for (int b=i+1; b<this->fsize; b++) {
            if (this->fvector [b] < this->fvector[menor]) {
                Swap(b, menor);
            }
        }
    }
    
}

/// doubleiza a troca dos elementos "a" e "small";
template<class T>
void TVecNum<T>::Swap (int a, int small)
{
    if( this->fvector == NULL)
        DebugStop();
    
    if (a<0 || a>= this->Size() || small < 0 || small >= this->Size()) {
        DebugStop();
    }
    T s = this->fvector[small];
    this->fvector[small]=this->fvector[a];
    this->fvector[a]=s;
}

///retorna o valor da forma do vetor criado na classe;
template<class T>
double TVecNum<T>::Norma() const
{
    if( this->fvector == NULL)
        DebugStop();
    
    double norma=0;
    
    for (int n=0; n<this->fsize; n++)
    {
        norma = norma + this->fvector[n] * this->fvector[n];
    }
    
    return sqrtf((double)norma);
}

/// zera os valores do vetor
template<class T>
void TVecNum<T>::Zero()
{
    for (int i=0; i<this->fsize; i++) {
        this->fvector[i] = 0;
    }
}



///cria-se um vetor vecSoma (variavel local), atribui-se a este vetor a soma do vetor vecFator(constante e passado por referencia) com o proprio vetor que chamou o operador +;
///eh retornado o vetor vecSoma;
template<class T>
TVecNum<T> TVecNum<T>::operator + (const TVecNum<T> &vecFator) const
{
    if (this->fsize != vecFator.Size()) {
        //cout << "Erro: soma de vetores de tamanhos diferentes.";
        DebugStop();
    }
    
    if (this->fvector == NULL) {
        DebugStop();
    }
    
    TVecNum<T> vecSoma(this->fsize);
    
    for (int i=0; i<this->fsize; i++) {
        vecSoma.fvector[i] = this->fvector[i] + vecFator.fvector[i];
    }
    return vecSoma;
}


///segue o mesmo raciocinio do operador +;
template<class T>
TVecNum<T> TVecNum<T>::operator - (const TVecNum<T> &vecFator) const
{
    if (this->fsize != vecFator.Size()) {
        //cout << "Erro: soma de vetores de tamanhos diferentes.";
        DebugStop();
    }
    
    if (this->fvector == NULL) {
        DebugStop();
    }
    
    TVecNum<T> vecDif(this->fsize);
    
    for (int i=0; i<this->fsize; i++) {
        vecDif.fvector[i] = this->fvector[i] - vecFator.fvector[i];
    }
    return vecDif;
}

///implementa produto interno
template<class T>
T TVecNum<T>::operator * (const TVecNum<T> &vecFator) const
{
    if (this->fsize != vecFator.Size()) {
        //cout << "Erro: soma de vetores de tamanhos diferentes.";
        DebugStop();
    }
    
    if (this->fvector == NULL) {
        DebugStop();
    }
    
    T  dotProduct = 0;
    
    for (int i=0; i<this->fsize; i++) {
        dotProduct += this->fvector[i] * vecFator.fvector[i];
    }
    return dotProduct;
}

template class TVecNum<int>;
template class TVecNum<double>;



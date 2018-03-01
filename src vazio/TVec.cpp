//  Vetor.cpp
//  tutorial_1
//
//  Created by Eduardo on 03/03/15.
//  Copyright (c) 2015 Eduardo. All rights reserved.
//

#include <iostream>
#include "TVec.h"
#include <stdio.h>
#include <math.h>
#include "tpanic.h"

#ifdef WIN32
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif


using std::cout;
using std::endl;
using std::cin;

//construtor padrao de TVec
template<class T>
TVec<T>::TVec()
{
    ///cout << __PRETTY_FUNCTION__ << " endereco " << (void *) this << std::endl;
    fvector = NULL;
    
    fsize = 0;
}

//fvector recebe o tamanho (fsize) em alocacao dinamica
//todos os elementos do vetor recebem zero
template<class T>
TVec<T>::TVec(int size)
{
    if (size < 0) {
        DebugStop();
    }
    this->fvector = new T[size];
    ///cout << __PRETTY_FUNCTION__ << " Endereco do vetor criado\t" << (void *) this << std::endl;
    fsize = size;
}

///construtor copia
template<class T>
TVec<T>::TVec(const TVec &cp)
{
    if( cp.fvector == NULL)
        DebugStop();
    
    fsize = cp.fsize;
    ///cout << __PRETTY_FUNCTION__ << " Endereco do vetor criado por copia\t" << (void *) this << std::endl;
    fvector = new  T [fsize];
    
    for (int i=0; i<fsize; i++) {
        fvector[i] = cp.fvector[i];
    }
}

///Destrutor dos vetores criados
template<class T>
TVec<T>::~TVec(void)
{
    ///cout << __PRETTY_FUNCTION__ << " Destruicao do vetor e liberacao de memoria\t" << (void *) this << endl;
    if(fvector != NULL){
        delete [] this->fvector;
    }
}


///Funcao do tipo Get; aponta para o a posicao e endereco do vetor;
template<class T>
T TVec<T>::GetVal(int i) const
{
    if (i >= fsize || i < 0) {
        //std::cout << "\nERRO - Acessando posicao inexistente" << std::endl;
        DebugStop();
    }
    return fvector[i];
}


///Funcao do tipo Put; permite o recebimento de valores e sua colocacao na respectiva parte do vetor;
template<class T>
void TVec<T>::PutVal(int i,  T val)
{
    if (i >= fsize || i < 0) {
        //std::cout << "\nERRO - Acessando posicao inexistente" << std::endl;
        DebugStop();
    }
    fvector[i] = val;
}




///cria um ponteiro "*copia" que recebe os valores do fvector, modifica o tamanho do vetor e "*copia" devolve os valores para o vetor de novo tamanho;
template<class T>
void TVec<T>::Resize (int nsize)
{
    if( nsize < 0 )
        DebugStop();
    //
    ///std::cout << "Expansao/Trucamento de vetor:";
    T *copia;
    copia = new  T [nsize];
    
    if(fvector!=NULL){
        for (int i=0; i<fsize; i++) {
            if (i<nsize)
            {
                copia[i]=fvector[i];
            }
        }
    }
    
    if(fvector != NULL){
        delete [] fvector;
    }
    fvector = copia;
    fsize = nsize;
}


///retorna o valor de fsize;
template<class T>
int	TVec<T>::Size() const
{
    return fsize;
}

///operador de selecao [], recebe a posicao "i" dentro do vetor e retorna o valor por referencia;
template<class T>
T &TVec<T>::operator [](int i)
{
    if(i<0 || i>=fsize)
    {	//cout << "Posicao inexistente no vetor" << endl;
        DebugStop();
    }
    
    return fvector[i];
}

///operador de selecao [], recebe a posicao "i" dentro do vetor e retorna o valor por referencia;
template<class T>
T &TVec<T>::operator [](int i) const
{
    if(i<0 || i>=fsize)
    {	//cout << "Posicao inexistente no vetor" << endl;
        DebugStop();
    }
    
    return fvector[i];
}


//operador =, recebe um vetor constante "vecCop" por referencia (nao ha alteracao nos valores do vecCop, passagem por referencia para nao ser necessario criar copia);
///o vetor que chamou o operador eh alterado por referencia (nao const), e recebe como retorno o proprio vetor que foi chamado, ja alterado;
template<class T>
TVec<T> &TVec<T>::operator = (const TVec<T> &vecCop)
{
    
    fsize = vecCop.fsize;
    if(fvector!=NULL)
    {
        delete [] fvector;
    }
    fvector = new  T [fsize];
    
    for (int i=0; i<fsize; i++) {
        fvector[i] = vecCop.fvector[i];
    }
    return *this;
}



///compara o vetor que chamou o operados com o vetor vec1 recebido por referencia;
///retorna se os vetores sao iguais ou nao;
template<class T>
bool TVec<T>::operator == (const TVec &vec1) const
{
    if(vec1.Size() != Size())
        return false;
    
    for (int i=0; i<fsize; i++){
        if(!( vec1.fvector[i] == fvector[i] ))
            return false;
    }
    
    return true;
}

//imprime vetor

template<class T>
void TVec<T>::Print(std::ostream &out) const
{
    for (int irow=0;irow<fsize; irow++) {
        out<<"O elemento ["<<irow<<"] é "<< GetVal(irow) <<std::endl;
    }
    
}



#include "tno.h"
class TElemento;

template class TVec<double>;
template class TVec<int>;
template class TVec<TNo>;
template class TVec<TElemento *>;


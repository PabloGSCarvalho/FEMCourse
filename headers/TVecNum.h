#ifndef TVECNUM_H
#define TVECNUM_H

/*
 *  TMatrix.h
 *
 *  Created by Eduardo Ferri on 3/6/15.
 *  Copyright 2015 LabMeC. All rights reserved.
 *
 */

#include <stdio.h>
#include <iostream>
#include "TVec.h"

template<class T>
class TVecNum : public TVec<T>
{
public:
  
  /// Construtor vazio
  TVecNum ();
  
  /// Construtor com tamanho definido do vetor
  TVecNum (int size);
  
  /// Construtor com tamanho definido do vetor
  TVecNum (int size, T value);
  
  /// Construtor de copia
  TVecNum (const TVecNum &cp);
  
  //Destrutor: zera o vetor
  ~TVecNum ();
  
  ///funcao que calcula norma
  double Norma() const;
  
  ///funcao que ordena os valores
  void OrdVec();
  
  ///funcao que troca os valores do vetor
  void Swap (int a, int small);
  
  ///funcao que expande o tamanho do vetor de fsize para nsize
  void Resize (int nsize);
    
  /// zera os valores do vetor
  void Zero();
  
  TVecNum operator + (const TVecNum &vecFator) const;
  
  TVecNum operator - (const TVecNum &vecFator) const;
  
  ///implementa produto interno entre dois vetores
  T operator * (const TVecNum &vecFator) const;
  
  
  
};

#endif
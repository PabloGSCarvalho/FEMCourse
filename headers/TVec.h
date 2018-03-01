#ifndef TVEC_H
#define TVEC_H

/*
 *  TVec.h
 *
 *  Created by Eduardo Ferri on 3/6/15.
 *  Copyright 2015 LabMeC. All rights reserved.
 *
 */

#include <stdio.h>
#include <iostream>

template<class T>
class TVec
{
public:
  
  /// Construtor vazio
  TVec ();
  
  /// Construtor com tamanho definido do vetor
  TVec (int size);
  
  /// Construtor de copia
  TVec (const TVec &cp);
  
  //Destrutor
  ~TVec ();
  
  ///funcao retorno os valores do vetor
  T GetVal(int i) const;
  
  ///funcao recebe os valores do vetor
  void PutVal(int i, T val);
  
  ///funcao que expande o tamanho do vetor de fsize para nsize
  void Resize (int nsize);
  
  
  ///funcao que retorna o tamanho do vetor
  int	Size() const;
  
  ///operador colchete para selecao de elemento quando lvalue
  T &operator [](int i);
  
  ///operador colchete para selecao de elemento quando rvalue
  T &operator [](int i) const;
  
  ///operador de atribuicao
  TVec &operator = (const TVec &vecCop);
  
  ///comparador de igualdade entre dois vetores
  bool operator == (const TVec &vec1) const;
  
  //imprime vetor
  void Print(std::ostream &out = std::cout) const;
  
  
protected:
  //ponteiro vetor
  T *fvector;
  //tamanho vetor
  int fsize;
};

#endif
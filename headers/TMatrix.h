#ifndef TMATRIX_H
#define TMATRIX_H

/*
 *  TMatrix.h
 *
 *  Created by Eduardo Ferri on 3/6/15.
 *  Copyright 2015 LabMeC. All rights reserved.
 *
 */

#include"TVecNum.h"
#include <iostream>
#include <stdio.h>

class TMatrix{
  
public:
  ///Construtor vazio;
  TMatrix ();
  
  
  ///Cosntrutor da matriz com row`s linhas e col`s colunas;
  TMatrix (int row, int col, double val = 0.);
  
  ///Destrutor;
  ~TMatrix ();
  
  ///Contrutor de cpoia. Recebe os valores de copia por referencia;
  TMatrix(const TMatrix &cp);
  
  ///Funcoes/Objetos
  
  ///Recebe o valor val e aloca o mesmo na posicao i j da matriz;
  void PutVal (int i,int j, double val);
  
  ///Retorna um valor double da matriz pertencenet a posicao i j da matriz;
  double GetVal (int i, int j) const;
  
  ///Retorna o numero de linhas que a matriz possui;
  int Rows(void) const;
  
  ///Retorna o numero de colunas que a matriz possui;
  int Cols(void) const;
  
  ///Edita o tamanho da matriz para um novo numero de row`s e col`s;
  void Resize(int row, int col);
  
  ///Zera os valores da matriz
  void Zero();
  
  ///Transpoe a matriz;
  void Transpose();
  
  ///Retorna um vetor do tamanho do numero de colunas da matriz referente a linha row;
  TVecNum<double> GetRow(int row) const;
  
  ///Retorna um vetor do tamanho do numero de linhas da matriz referente a coluna col;
  TVecNum<double> GetCol(int col) const;
  
  /// Decompoe a matriz
  void LU_Decomposition();

  ///Resolve o sistema linear
  void Solve_LU(TMatrix &rhs);
  
  ///Imprime a matriz
  void Print(std::ostream &out = std::cout)
  {
    for (int irow=0;irow<Rows(); irow++) {
      //out << "Linha " << irow << " : ";
      for (int icol=0;icol<Cols(); icol++) {
        out<<GetVal(irow, icol) << " ";
      }
      out << std::endl;
    }
  }
  
  //Imprime a matriz
  void Print(const std::string &txt, std::ostream &out = std::cout)
  {
    out << txt.c_str() << std::endl;
    Print(out);
  }
    
//   static bool IsZero( double val);
    
    
  ///Operadores
  
  ///Chama o valor existente na posicao i,j da matriz, sendo possivel, retornar o valor existente nesta posicao ou atribuir um valor double a esta posicao;
  double &operator ()(int i, int j);
  
  ///A matriz que chamou recebe por referencia a matriz copia (similar ao construtor de copia);
  TMatrix &operator =(const TMatrix &copia);
  
  ///Soma os elementos de duas matrizes de mesmo tamanho, retornando uma matriz com este valor;
  TMatrix operator +(TMatrix &fator);
  
  ///Subtrai os elementos de duas matrizes de mesmo tamanho, retornando uma matriz com este valor;
  TMatrix operator -(TMatrix &fator);
  
  ///Multiplica todos os elementos de uma matriz por um double m;
  TMatrix &operator *(double m);
  
  ///Multiplica duas matrizes respeitando as regras da algebra linear;
  TMatrix operator *(TMatrix &fator);
  
  ///Multiplica uma matriz por um vetor vec e retorna um valor de mesmo tamanho vec (respeita as regras da algebra linear);
  TVecNum<double> operator *(TVecNum<double> &vec);
  
  ///Compara duas matrizes e retorna true ou false;
  bool operator == (TMatrix &mat2);
  
private:
  /// ponteiro de uma dimensao que aponta para os elementos constituintes da matriz;
  double *fmatrix;
  
  ///numero de linhas da matriz
  int frow;
  int fcol;
  int fdecompose;
};

std::ostream &operator<<(std::ostream &out, const TMatrix &mat);

#endif
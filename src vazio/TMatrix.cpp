/*
 *  TMatrix.cpp
 *  Tutorial_1
 *
 *  Created by Eduardo Ferri on 3/6/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "TMatrix.h"
#include "tpanic.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

using std::cout;
using std::endl;
using std::cin;

TMatrix::TMatrix()
{
  
  fmatrix= NULL;
  frow =0;
  fcol=0;
  fdecompose=0;
  
}

TMatrix::TMatrix(int row, int col, double val)
{
  if (row < 0 || col < 0) {
    DebugStop();
  }
  frow = row;
  fcol = col;
  fdecompose=0;
  fmatrix = new double [frow*fcol];
  for (int i=0; i<frow*fcol; i++) {
    fmatrix[i] = val;
  }
}

TMatrix::~TMatrix(void)
{
  if(fmatrix != NULL){
    for (int i=0; i<frow*fcol; i++) {
      fmatrix[i] = 0.;
    }
    delete [] this->fmatrix;
  }
}

TMatrix::TMatrix(const TMatrix &cp)
{
  
  frow= cp.frow;
  fcol= cp.fcol;
  fdecompose=cp.fdecompose;
  
  fmatrix = new double [frow*fcol];
  
  for (int i=0; i<frow*fcol; i++) {
    fmatrix[i] = cp.fmatrix[i];
  }
}

double TMatrix::GetVal (const int i,const int j) const
{
  if (i<0 || i>=frow || j<0 || j>=fcol) {
    DebugStop();
  }
  return fmatrix[frow*j + i];
}

int TMatrix::Rows(void) const
{
  return frow;
}
int TMatrix::Cols(void) const
{
  return fcol;
}

void TMatrix::PutVal (int i, int j, double val)
{
  if (i<0 || i>=frow || j<0 || j>=fcol) {
    DebugStop();
  }
  fmatrix[frow*j + i] = val;
}

void TMatrix::Resize(int row, int col)
{
  if ( row <= 0 || col <= 0 ||fdecompose==1) {
    DebugStop();
  }
  
  if (fmatrix == NULL) {
    fmatrix = new double[row*col];
    for (int i=0; i<row*col; i++) {
      fmatrix[i] = 0.;
    }
    frow = row;
    fcol = col;
  }
  else{
    TMatrix copy(*this);
    
    
    delete [] fmatrix;
    
    fmatrix = new double[row*col];
    frow=row;
    fcol=col;
    
    for (int i=0; i<row; i++) {
      for (int j=0; j<col; j++){
        if (i>=copy.Rows())
        {
          PutVal(i, j, 0);
        }
        else if (j>=copy.Cols())
        {
          PutVal(i, j, 0);
        }
        else
        {
          PutVal(i, j, copy(i,j));
        }
      }
    }
  }
}


void TMatrix::Transpose()
{
  TMatrix tr(fcol,frow);
  
  if (fmatrix==NULL||fdecompose==1) {
    DebugStop();
  }
  
  for (int i=0; i<frow; i++) {
    for (int j=0; j<fcol; j++) {
      tr.PutVal(j, i, GetVal(i, j));
    }
  }
  fcol=tr.fcol;
  frow=tr.frow;
  
  for (int i=0; i<frow*fcol; i++) {
    this->fmatrix[i]=tr.fmatrix[i];
  }
}

///Zera os valores da matriz
void TMatrix::Zero()
{
  int sz = frow*fcol;
  for (int i=0; i<sz; i++) {
    fmatrix[i] = 0.;
  }
  fdecompose=0;
  
}


double &TMatrix::operator ()(int i, int j)
{
  if ( i< 0 || i >= frow ) {
    DebugStop();
  }
  if( j<0 || j >= fcol ) {
    DebugStop();
  }
  return fmatrix[frow*j + i];
}


TMatrix &TMatrix::operator=(const TMatrix &copia)
{
  
  if(&copia!=this){
    if(frow!=copia.frow||fcol!=copia.fcol){
      delete []fmatrix;
      frow=copia.frow;
      fcol=copia.fcol;
      fdecompose=copia.fdecompose;
      fmatrix=new double[frow*fcol];
      
    }
    
    for (int i=0; i<frow*fcol; i++) {
      fmatrix[i]=copia.fmatrix[i];
    }
    
  }
  
  return *this;
  
}

TMatrix TMatrix::operator +(TMatrix &fator)
{
  
  
  if(frow!=fator.Rows() || fcol!=fator.Cols()||fator.fdecompose==1||fdecompose==1)
  {
    cout << "Nao eh possivel realizar a soma: matrizes de tamanho diferentes." << endl;
    
  }
  
  TMatrix soma(frow,fcol);
  
  for (int i=0; i<frow; i++) {
    for(int j=0;j<fcol;j++){
      soma.PutVal(i, j, (GetVal(i, j)+fator.GetVal(i, j)));
      ///soma.fmatrix[i]=this->fmatrix+fator.fmatrix[i];
    }
  }
  return soma;
}

TMatrix TMatrix::operator -(TMatrix &fator)
{
  if(frow!=fator.Rows() || fcol!=fator.Cols()||fator.fdecompose==1||fdecompose==1)
  {
    cout << "Nao eh possivel realizar a soma: matrizes de tamanho diferentes." << endl;
    
  }
  
  TMatrix dif(frow,fcol);
  
  for (int i=0; i<frow; i++) {
    for(int j=0;j<fcol;j++){
      dif.PutVal(i, j, (GetVal(i, j)-fator.GetVal(i, j)));
    }
  }
  return dif;
}

TMatrix &TMatrix::operator *(double m)
{
  if (fdecompose==1) {
    DebugStop();
    
  }
  
  for (int i=0; i<frow*fcol; i++) {
    fmatrix[i]=m*fmatrix[i];
  }
  
  return *this;
}

///PRODUTO DE MATRIZES
TMatrix TMatrix::operator *(TMatrix &fator)
{
  if (fdecompose==1) {
    DebugStop();
		}
  
  if (fcol!=fator.Rows()) {
    DebugStop();
    //cout << "O numero de colunas da matriz A eh diferente do numero de linhas da matriz matriz vazia. Nao eh possivel realizar a multiplicacao" <<endl;
    
  }
  
  TMatrix mult(frow,fator.Cols());
  for (int i=0; i<mult.Rows(); i++) {
    for (int j=0; j<mult.Cols(); j++) {
      for (int k=0; k<Cols(); k++) {
        mult(i, j) += GetVal(i, k)*fator(k, j);
      }
    }
  }
  
  return mult;
}


///PRODUTO MATRIZ VETOR
TVecNum<double> TMatrix::operator *(TVecNum<double> &vec)
{
  if (fcol!=vec.Size()){
    DebugStop();//cout << "Numero de colunas da matriz diferente do tamanho do vetor." <<endl;
    
  }
  
  TVecNum<double> mult(frow);
  
  for (int j=0; j<frow; j++) {
    for (int i=0; i<fcol; i++) {
      mult[j] += GetVal(j,i)*vec[i];
    }
  }
  
  return mult;
}

///OPERADOR ==
bool TMatrix::operator == (TMatrix &mat2)
{
  if (mat2.fmatrix == NULL) {
    DebugStop();
  }
  if(mat2.Rows() != Rows() || mat2.Cols() != Cols())
    return false;
  
  for (int i=0; i<frow; i++){
    for (int j=0; j<fcol; j++) {
      if(mat2.GetVal(i, j) != GetVal(i, j))
        return false;
    }
  }
  
  return true;
}

//static bool IsZero( double val)
//{
//  return fabs(val) < 1.e-9;
//}

//Decmopoe a matriz

void TMatrix::LU_Decomposition()
{
  
	if (fdecompose == 1){
		cout << "Decomposicao ja foi realizada" << endl;
		DebugStop();

	}

	int s;
	double aux = 0.;

	if (Rows() != Cols()) {
		DebugStop();
	}


	if (Rows()<Cols()) {
		s = Rows();
	}
	else{
		s = Cols();
	}

	for (int k = 0; k<s; k++)
	{
		if (GetVal(k, k) == 0) {
			DebugStop();

		}
		else{

			for (int irow = k + 1; irow<Rows(); irow++) {
				aux = GetVal(irow, k) / GetVal(k, k);
				PutVal(irow, k, aux);

				for (int icol = k + 1; icol<Cols(); icol++) {
					PutVal(irow, icol, GetVal(irow, icol) - GetVal(k, icol)*aux);
				}

			}


		}

	}

	fdecompose = 1;
  
}


//Resolve o sistema
void TMatrix::Solve_LU(TMatrix &rhs)
{
  
	if (fmatrix == NULL) {
		DebugStop();
	}

	LU_Decomposition();

	if (rhs.Rows() != Cols()) {
		DebugStop();
	}

	//Resolucao do sistema

	for (int irow = 0; irow<Rows(); irow++) {

		for (int icol = 0; icol<rhs.Cols(); icol++) {

			double aux2 = 0.;

			for (int k = 0; k<irow; k++) {
				aux2 += rhs.GetVal(k, icol)*GetVal(irow, k);
			}

			rhs.PutVal(irow, icol, rhs.GetVal(irow, icol) - aux2);


		}


	}


	for (int irow = (Rows() - 1); irow >= 0; irow--) {

		for (int icol = 0; icol<rhs.Cols(); icol++) {

			double aux1 = 0.;

			for (int k = (irow + 1); k<Cols(); k++) {

				aux1 += rhs.GetVal(k, icol)*GetVal(irow, k);
			}


			rhs.PutVal(irow, icol, (rhs.GetVal(irow, icol) - aux1) / GetVal(irow, irow));


		}
	}


  
}

TVecNum<double> TMatrix::GetRow(int row) const {
	
	if (row<0 || row >= Rows()){
		DebugStop();
	}

	TVecNum<double> VparcC(Cols());
	for (int icol = 0; icol<Cols(); icol++){

		VparcC[icol] = GetVal(row, icol);
	}


	return VparcC;

}

TVecNum<double> TMatrix::GetCol(int col) const {
	
	if (col<0 || col >= Cols()){
		DebugStop();
	}

	TVecNum<double> VparcR(Rows());
	for (int irow = 0; irow<Rows(); irow++){

		VparcR[irow] = GetVal(irow, col);
	}


	return VparcR;
}

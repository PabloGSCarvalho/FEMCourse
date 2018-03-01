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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <cstdlib>
#include <math.h>
#include "tno.h"
#include "telemento.h"
#include "telemento0d.h"
#include "telemento1d.h"
#include "tmalha.h"
#include "tanalysis.h"
#include "tmaterial1d.h"
#include "tmaterialbc.h"

//using namespace std;

double Integrate(int ord, double xmin, double xmax,double (*g)(double x));

double func(double x);

void ReadMesh(TMalha &malha, std::string &filename);

/*
int main(int argc, char *argv[])
{
  std::cout << "Hello, world!" << std::endl;
  
  std::cout << "Testando classe nó" << std::endl;
  TNo::main();
  std::cout << "Fim do teste da classe nó" << std::endl;

  std::cout << "Testando classe elemento" << std::endl;
  TElemento::main();
  std::cout << "Fim do teste da classe elemento" << std::endl;
  
  std::cout << "Testando uma malha uni dimensional" << std::endl;
  TMalha::main();
  std::cout << "Fim do teste da class malha" << std::endl;
  
  double result;
  int ord;
  for (ord = 1; ord < 10; ord++)
  {
    result = Integrate(ord,0.,3.,func);
    std::cout << "O resultado da integracao " << result << std::endl;
  }

  if(argc < 2)
  {
      std::cout << "Numero de argumentos menor que dois";
        return -1;
  }
  TMalha NumeroUm;
  std::string filename(argv[1]);
  ReadMesh(NumeroUm,filename);
  TAnalysis analysis(&NumeroUm);
  analysis.Run();
    
  return EXIT_SUCCESS;
}

double Integrate(int ord, double xmin, double xmax,double (*g)(double x))
{
  // escala devido ao tamanho do intervalo
  double mult = (xmax-xmin)/2.;
  // objecto tipo regra de integracao
  TPZInt1d rule(ord);
  // numero de pontos de integracao
  int np = rule.NPoints();
  std::cout << "Numero de pontos " << np << std::endl;
  // estrutura de dados para receber a posicao do ponto e o peso
  TPZVec<double> point(1);
  double w, result = 0.;
  int p;
  // loop sobre os pontos de integracao
  for(p=0; p<np; p++)
  {
    // pegue a posicao do ponto e o peso
    rule.Point(p,point,w);
    // calcula a coordenada x correspondente
    double x = xmin + (point[0]+1.)*(xmax-xmin)/2.;
    // calcula o valor da funcao, integra
    result += g(x)*w*mult;
    
  }
  return result;
  
}

double func(double x)
{
  return sin(x);
}

void ReadMesh(TMalha &malha, std::string &filename)
{
   int nmat,nbc,npoint;
   std::ifstream input(filename.c_str());
   input >> nmat >> nbc >> npoint;
   int im;
   for(im=0; im<nmat; im++)
   {
   	int id;
	double k,c,b,f;
	input >> id >> k >> c >> b >> f;
	TMaterial1d *mat1d = new TMaterial1d(id,k,c,b,f);
	malha.insertMaterial(mat1d);
   }
   for(im=0; im<nbc; im++)
   {
   	int id, type;
	double contrstiff, contrrhs;
	input >> id >> type >> contrstiff >> contrrhs;
	TMaterialBC *matbc = new TMaterialBC(id,type,contrstiff,contrrhs);
	malha.insertMaterial(matbc);
   }
   for(im=0; im < npoint; im++)
   {
   	int mat;
	int node, order = 0;
   	input >> mat >> node;
	std::vector<int> nodes(1,node);
	TElemento0d *elem = new TElemento0d(mat,order,nodes);
	malha.insertElement(elem);
   }
   double x0;
   int matleft, matright;
   input >> x0 >> matleft >> matright;
   std::vector<double> co(1,x0);
   TNo no(co);
   malha.insertNode(no);
   int count = 0;
   while(input)
   {
   	double x1;
	int nelem, mat, morder;
	input >> x1 >> nelem >> mat >> morder;
	if(!input) break;
	int numnewnodes = nelem*morder;
	int firstnode = malha.getNodeVec().size();
	malha.getNodeVec().resize(firstnode+numnewnodes);
	int in;
	for(in=1; in<=numnewnodes; in++)
	{
		co[0] = x0+in*(x1-x0)/numnewnodes;
		TNo no(co);
		malha.getNode(firstnode+in-1) = no;
	}
	int el;
	count = firstnode-1;
	for(el=0; el<nelem; el++)
	{
		std::vector<int> nodes(morder+1);
		int in;
		for(in=0; in<=morder; in++)
		{
			nodes[in] = count++;
		}
		TElemento1d *elem = new TElemento1d(mat,morder,nodes);
		malha.insertElement(elem);
		count--;
	}
   }
   int node, order = 0;
   std::vector<int> nodes(1,0);
   TElemento0d *elem = new TElemento0d(matleft,order,nodes);
   malha.insertElement(elem);
   nodes[0] = count;
   elem = new TElemento0d(matright,order,nodes);
   malha.insertElement(elem);
}
  
  
  }
 */

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


#include "telemento1d.h"

#include "tno.h"
#include "tmalha.h"
#include "tmaterial.h"
#include "TMatrix.h"
#include "TIntRule1d.h"
#include "tpanic.h"
#include <math.h>
//#include "pzquad.h"

TElemento1d::TElemento1d(int matid, int order, TVec< int > &nodes): TElemento(matid, order, nodes)
{
}


TElemento1d::TElemento1d(): TElemento()
{
}


TElemento1d::~TElemento1d()
{
}


MElementType TElemento1d::getElType()
{
    return ELinear;
}

void TElemento1d::CalcStiff(TMalha &malha, TMatrix &stiff, TMatrix &rhs)
{
    
    
    //identificar o objeto material
    
    int matId =this->fMaterialId;
    TMaterial *mat=malha.getMaterial(matId);
    
    //inicializar as matrizes
    
    stiff.Resize(fPorder+1, fPorder+1);
    rhs.Resize(fPorder+1, 1);
    
    //Criar regra de integracao
    
    TIntRule1d iRule1d(fPorder+fPorder);
    int np=iRule1d.NPoints();
    
    TMatrix dphi(1,fPorder+1);
    TVec<double> phi(fPorder+1);
    
    
    for (int i=0; i<np; i++) {
        double weight=0.0;
        TVec<double> co(1);
        iRule1d.Point(i,co,weight);
        Shape(co,phi,dphi);
        
        TMatrix jacobian(1,1);
        TMatrix jacinv(1,1);
        double detjac=0.0;
        
        Jacobian(co,jacobian,jacinv,detjac,malha);
        
        weight=weight*fabs(detjac);
        for (int i = 0; i<phi.Size(); i++) {
            //phi[i]=phi[i]*detjac;
            
            dphi(0,i)=dphi(0,i)*jacinv(0,0);
        }
        
        
        mat->Contribute(weight,phi,dphi,stiff,rhs);
       // stiff.Print("Pos Contribute");
    }
    
    
    
}


/**
 * Calculo do jacobiano
 * @point : ponto em qual calculamos o jacobiano do mapeamento
 * @jacobian : matriz jacobiana (a dimensao depende da dimensao do elemento
 * @jacinv : inverso da matriz jacobiana
 * @malha : objeto malha necessaria para relacionar os indices dos nos com os nos reais
 */

void TElemento1d::Jacobian(TVec<double> &point, TMatrix &jacobian, TMatrix &jacinv, double &detjac, TMalha &malha)
{
  
	TVec<TNo> nos(fNodes.Size());
	for (int ino = 0; ino<nos.Size(); ino++) {
		nos[ino] = malha.getNode(fNodes[ino]);
	}


	TMatrix Coord(2, fNodes.Size());

	for (int i = 0; i<fNodes.Size(); i++) {
		Coord(0, i) = nos[i].Co(0);
		Coord(1, i) = nos[i].Co(1);
	}

	TVec<double> phixi(fNodes.Size());
	TMatrix dphixi(1, fNodes.Size());

	jacobian.Resize(1, 1);
	jacobian.Zero();

	jacinv.Resize(1, 1);
	jacinv.Zero();

	this->Shape(point, phixi, dphixi);
	double dxdxi = 0., dydxi = 0.;
	for (int xi = 0; xi<fNodes.Size(); xi++) {
		dxdxi += Coord(0, xi)*dphixi(0, xi);
		dydxi += Coord(1, xi)*dphixi(0, xi);

	}

	jacobian(0, 0) = sqrt(dxdxi*dxdxi + dydxi*dydxi)/2.;

	jacobian.Resize(1, 1);

	jacinv(0, 0) = 1./jacobian(0,0);

	detjac = jacobian(0, 0);

	//

	/*
    TNo n01=malha.getNode(fNodes[0]);
    TNo n02=malha.getNode(fNodes[fPorder]);
    
    TVecNum<double> Coord0(2);
    TVecNum<double> Coord1(2);
    
    Coord0[0]=n01.Co(0);
    Coord0[1]=n01.Co(1);
    Coord0[0]=n02.Co(0);
    Coord0[1]=n02.Co(1);
    
    double hel = sqrt(
                      pow(n01.Co(0) - n02.Co(0), 2.0) +
                      pow(n01.Co(1) - n02.Co(1), 2.0));
    
    jacobian.Resize(1, 1);
    jacinv.Resize(1, 1);
    jacobian(0,0)=hel/2.0;
    jacinv(0,0)=2.0/hel;
    detjac=hel/2.0;
	*/

    
}

/**
 * Calcula os valores das funcoes de forma e suas derivadas
 * @point ponto onde calcular as funcoes de forma
 * @phi valores das funcoes de forma
 * @dphi valores das derivadas das funcoes de forma */


void TElemento1d::Shape(TVec<double> &point, TVec<double> &phi, TMatrix &dphi)
{
    
    for (int i=0; i<fPorder+1; i++) {
        phi[i]=1.;
        dphi(0,i)=0;
    }
    
    for (int i=0; i<fPorder+1; i++) {
        double epsi=-1.+i*2./fPorder;
        
        for (int j=0; j<fPorder+1; j++) {
            
            TMatrix axdphi(1,fPorder+1);
            if (i!=j) {
                double epsj=-1.+j*2./fPorder;
            
                phi[i]*=(point[0]-epsj)/(epsi-epsj);
                
                axdphi(0,i)=1/(epsi-epsj);
                
                for (int k=0; k<fPorder+1; k++) {
                    if (k!=i&&k!=j) {
                        epsj=-1.+k*2./fPorder;
                        axdphi(0,i)*=(point[0]-epsj)/(epsi-epsj);
                    }

                }
            
                dphi(0,i)+=axdphi(0,i);
            
            }
        
        
        }
    
    }
    
}

/**
 * Calcula o erro do elemento
 * @param exact funcao que calcula o valor exato
 * @param energy [out] erro na norma de energia
 * @param l2 [out] erro na norma l2
 */
void TElemento1d::Error(TMatrix &solution, TMalha &malha, void (exact)(TVec<double> &,double &, TVec<double> &), double &energy, double &l2)
{

    int matId =this->fMaterialId;
    TMaterial *mat=malha.getMaterial(matId);
    
    //calcluar as coordenadas , limites dos elementos
    TNo n01=malha.getNode(fNodes[0]);
    TNo n02=malha.getNode(fNodes[fPorder]);
    int nnodes=fNodes.Size();
    
    TVecNum<double> xA(1);
    TVecNum<double> xB(1);
    
    xA[0]=n01.Co(0);
    xB[0]=n02.Co(0);
    
    
    //calcular o número de pontos de integracao
    TIntRule1d iRule1d(19);
    int np=iRule1d.NPoints();
    energy=0.;
    l2=0.;

    for (int ip=0; ip<np; ip++) {
        //calcualndo o psi-ponto de integracao
        double weight=0.0;
        TVec<double> co(1);
        
        iRule1d.Point(ip, co, weight);
        //calculando o vetor point das cordenadas em funcao de psi (volta)
        
        TVec<double> xPoint(1);
        xPoint[0]=xA[0]+(xB[0]-xA[0])*(co[0]+1.0)/2.0;
        
        //calculo do jacobiano
        
        TMatrix jacobian(1,1);
        TMatrix jacinv(1,1);
        double detjac=0.0;
        
        Jacobian(co,jacobian,jacinv,detjac,malha);
        
        //calculo das funcoes de forma
        
        TMatrix dphi(1,fPorder+1);
        TVec<double> phi(fPorder+1);
        Shape(co,phi,dphi);
        
        weight=weight*fabs(detjac);
        for (int i = 0; i<phi.Size(); i++) {
            
            dphi(0,i)=dphi(0,i)*jacinv(0,0);
        }
        
        //calculo dos vetor de solucao e alphas
        
        double uh=0.;
        TVec<double> duh(1);
        duh[0]=0.;
        
        for (int i=0; i<nnodes; i++) {
            uh+=solution(fNodes[i],0)*phi[i];
            duh[0]+=solution(fNodes[i],0)*dphi(0,i);
        }
        
        mat->ContributeErrorSquare(xPoint, weight, uh, duh,
                                   exact, energy,l2);
        
    }
    
//    energy=sqrt(energy);
//    l2=sqrt(l2);
    
}

void TElemento1d::uhe(TMatrix &solution, TMalha &malha, TMatrix &uhe, TMatrix &duhedx){
    
    int matId =this->fMaterialId;
    TMaterial *mat=malha.getMaterial(matId);
    
    //calcluar as coordenadas , limites dos elementos
    TNo n01=malha.getNode(fNodes[0]);
    TNo n02=malha.getNode(fNodes[fPorder]);
    int nnodes=fNodes.Size();
    
    TVecNum<double> xA(1);
    TVecNum<double> xB(1);
    
    xA[0]=n01.Co(0);
    xB[0]=n02.Co(0);
    
    
    //calcular o número de pontos de integracao
//    TIntRule1d iRule1d(fPorder+fPorder);
//    int np=iRule1d.NPoints();

    int np = 3;
    uhe.Resize(np,2);
    duhedx.Resize(np, 2);
    uhe.Zero();
    duhedx.Zero();
    
    
    TVec<double> xi(3);
    
    xi[0]=-1.0;
    xi[1]= 0.0;
    xi[2]= 1.0;
    
    for (int ip=0; ip<np; ip++) {
        TVec<double> par(1);
        par[0]=xi[ip];
        
        TVec<double> xPoint(1);
        xPoint[0]=xA[0]+(xB[0]-xA[0])*(par[0]+1)/2;
        
        //calculo do jacobiano
        
        TMatrix jacobian(1,1);
        TMatrix jacinv(1,1);
        double detjac=0.0;
        
        Jacobian(par,jacobian,jacinv,detjac,malha);
        
        //calculo das funcoes de forma
        
        TMatrix dphi(1,fPorder+1);
        TVec<double> phi(fPorder+1);
        Shape(par,phi,dphi);
        
        for (int i = 0; i<phi.Size(); i++) {
            
            dphi(0,i)=dphi(0,i)*jacinv(0,0); // transforma a derivada no espaco parametrico para o dominio xyz. 
        }
        
        //calculo dos vetor de solucao gradiente e alphas
        
        uhe(ip,0)=xPoint[0];
        duhedx(ip,0)=xPoint[0];
        
        for (int i=0; i<nnodes; i++) {
            uhe(ip,1)+=solution(fNodes[i],0)*phi[i];
            duhedx(ip,1)+=solution(fNodes[i],0)*dphi(0,i);
        }
        
    }
    
    return;
}

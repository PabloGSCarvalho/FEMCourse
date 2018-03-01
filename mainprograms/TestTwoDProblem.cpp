//
//  TestOneDProblem.cpp MODIFICADO DO ORIGINAL
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
 Os testes foram preparados com um proposito educacional,
 recomenda-se que o aluno entenda a funcionalidade de cada
 teste e posteriormente use com seu cÛdigo caso a caso
 */
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>
#include "TMatrix.h"
#include "tmalha.h"
#include "telemento1d.h"
#include "telemento0d.h"
#include "telementoQuad.h"
#include "telementoTriangle.h"
#include "tmaterial1d.h"
#include "tmaterial2d.h"
#include "tmaterialbc.h"
#include "tanalysis.h"
#include "TIntRule1d.h"

#ifdef WIN32
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif

using std::cout;
using std::endl;
using std::cin;

void CreateTestMesh(TMalha &mesh, int order);
int TestNodes(TMalha *mesh);
int TestElement(TMalha *mesh);
int TestMesh(TMalha *mesh);

void CreateTestMesh(TMalha &mesh, int order);
void CreateTestMeshTR(TMalha &mesh, int order);
void TestOneDProblem(TMalha *mesh);

void exact(TVec<double> &point,double &val, TVec<double> &deriv);

int main ()
{
    
    
    TMalha mesh;
    int order = 2;
    //double h=1.; //para dar el tamanho del elemento
    CreateTestMesh(mesh,order);
	
    mesh.Print();
    
	TAnalysis Analisis(&mesh);
	Analisis.Run();
    
    double energy=0.0,l2=0.0;
    
    Analisis.Error(exact, energy, l2);
    
    cout<<endl;
    cout<<energy<<endl;
    cout<<l2<<endl;

    getchar();
    return 0;
}

void exact(TVec<double> &point,double &val, TVec<double> &deriv){
    
    
    long double Pi=4*atan(1);
    int mmax=200;
    int nmax=200;
    
//    val=0.;
//    deriv[0]=0.;
//    deriv[1]=0.;
//    
//    val=-point[0]*point[0]+point[0];
//    deriv[0]=-2.*point[0]+1.;
//    deriv[1]=0.;
    
    for(int m=1;m<=mmax;m++){
        for(int n=1;n<=nmax;n++){
            val+=(4*(-1 + pow(-1,m))*(-1 + pow(-1,n))*sin(m*Pi*point[0])*sin(n*Pi*point[1]))/(m*n*(pow(m,2) + pow(n,2))*pow(Pi,4));
            
            deriv[0]+=(4*(-1 + pow(-1,m))*(-1 + pow(-1,n))*cos(m*Pi*point[0])*sin(n*Pi*point[1]))/(n*(pow(m,2) + pow(n,2))*pow(Pi,3));
            
            deriv[1]+=(4*(-1 + pow(-1,m))*(-1 + pow(-1,n))*cos(n*Pi*point[1])*sin(m*Pi*point[0]))/(m*(pow(m,2) + pow(n,2))*pow(Pi,3));
            
        }
    }
    
    
    
    
}

void CreateTestMesh(TMalha &mesh, int order)
{
    double xmin = 0.;
    double xmax = 1.;
    double ymin = 0.;
    double ymax = 1.;
    
    //int nelemx = (xmax-xmin)/h;
    //int nelemy = (ymax - ymin) / h;
    
	int nelemx = 32;
	int nelemy = nelemx;   //mismo numero de elementos
    
    int nelem = nelemx*nelemy;
    
    
    
    int elorder = order;
    int nnodesx = nelemx*elorder + 1;  //nodes en x
    int nnodesy = nelemy*elorder + 1;  //mismo numero de nodes
	int nnodes = nnodesx*nnodesy;
    
	int matid = 1;
    double K = 1., B = 0., C[2] = { 0, 0 }, F(1.);  //elemento2d
    
	int leftmatid = -1;
    int rightmatid = -2;
	int downmatid = -3;
    int upmatid = -4;
    
    
    
    TVec<TNo> &nodevec = mesh.getNodeVec();
    nodevec.Resize(nnodes);
    
	TVecNum<double> coordx(nnodesx);
	TVecNum<double> coordy(nnodesy);
    
    for (int i=0; i<nnodesx; i++) {
        coordx[i] = xmin + (xmax-xmin)*i*1./(nnodesx-1);
	}
    
    for (int j=0; j<nnodesy; j++) {
        coordy[j] = ymin + (ymax-ymin)*j*1./(nnodesy - 1);
    }
    
    TVecNum<double> coord(2);
    
    for (int i = 0; i < nnodesx; i++)
    {
        for (int j = 0; j < nnodesy; j++)
        {
			coord[0] = coordx[j];
			coord[1] = coordy[i];
			nodevec[nnodesy*i+j].setData(coord);
        }
    }
    
    
    
    
    TVec<TElemento *> &elvec = mesh.getElementVec(); //inicializacao do objeto
    //elvec.Resize(nelem+4); //nelm elementos de interior + 2 elementos de BC;
    //elvec.Resize(nelem); //por lo pronto sin elementos de contorno;
	int nelem1D=nelemx*2+nelemy*2;
	elvec.Resize(nelem+nelem1D);
	
	
	TVec<int> nodes((elorder+1)*(elorder+1)); //node index no interior do elemento;
	
	for (int ely=0; ely<nelemy; ely++) {
		for (int elx = 0; elx<nelemx; elx++)
		{
            
            switch (elorder){
                case 1:
                    
                    nodes[0]=0;					//para el elemento 0 orden 1
                    nodes[1]=1;
                    nodes[2]=nnodesx+1;
                    nodes[3]=nnodesx;
                    
                    for (int i=0; i<((elorder+1)*(elorder+1)); i++)
                    {
                        nodes[i]= nodes[i]+elx+ely*nnodesx;
                        
                    }
                    
                    break;
                case 2:
                    
                    nodes[0]=0;					//para el elemento 0 orden 1
                    nodes[1]=2;
                    nodes[2]=(nnodes-1)/nelemx;
                    nodes[3]=((nnodes-1)/nelemx)-2;
                    nodes[4]=1;
                    nodes[5]=nnodesx+2;
                    nodes[6]=nnodesx+nnodesy+1;
                    nodes[7]=nnodesx;
                    nodes[8]=nnodesx+1;
                    
                    
                    for (int i=0; i<((elorder+1)*(elorder+1)); i++)
                    {
                        nodes[i]= nodes[i]+2*elx+2*ely*nnodesx;
                    }
					
                    break;
                case 3:
                    
                    nodes[0]=0;					//para el elemento 0 orden 1
                    nodes[1]=3;
                    nodes[2]=(nnodes-1)/nelemx;
                    nodes[3]=((nnodes-1)/nelemx)-3;
                    nodes[4]=1;
                    nodes[5]=2;
                    nodes[6]=nnodesx+3;
                    nodes[7]=nnodesx+nnodesy+3;
                    nodes[8]=(nnodesx+nnodesy+3)+nelemx*3;
                    nodes[9]=(nnodesx+nnodesy+3)+nelemx*3-1;					//para el elemento 0 orden 1
                    nodes[10]=nnodesx+nnodesy;
                    nodes[11]=nnodesx;
                    nodes[12]=nnodesx+1;
                    nodes[13]=nnodesx+2;
                    nodes[14]=nnodesx+nnodesy+2;
                    nodes[15]=nnodesx+nnodesy+1;
                    
                    
                    for (int i=0; i<((elorder+1)*(elorder+1)); i++) 
                    {
                        nodes[i]= nodes[i]+3*elx+3*ely*nnodesx;
                        
                    }
                    
                    break;
                }
			elvec[nelemy*ely+elx] = new TElementoQuad(matid,elorder,nodes);
		}
        
    }
    
	TMaterial *mat = new TMaterial2d(matid,K,C,B,F);
    mesh.insertMaterial(mat);
    
	
	/*TMatrix stiff, rhs;
     mesh.getElement(0)->CalcStiff(mesh,stiff,rhs);
     std::cout<<"matriz rigidez elemeto 0"<<std::endl;
     stiff.Print();
     */
	
	
	//condicoes de contorno
    
	TVecNum<double> nodesxdown(nnodesx);
    TVecNum<double> nodesxup(nnodesx);
	TVecNum<double> nodesyleft(nnodesy);
	TVecNum<double> nodesyright(nnodesy);
    
    for (int i=0; i<nnodesx; i++)
    {
        nodesxdown[i]= i;
        nodesyleft[i]=i*nnodesx;
        nodesxup[i]= i+nnodesx*(nnodesy-1);
        nodesyright[i]=i*nnodesx+(nnodesx-1);
    }
    
    
	
	
    TVec<int> nodesbc(elorder+1);
    
    for (int elbc=0; elbc<nelemx; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesxdown [elbc*(elorder)+i];
        }
        elvec[nelem+elbc] = new TElemento1d(downmatid,elorder,nodesbc);
    }
    
    for (int elbc=0; elbc<nelemx; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesxup [elbc*(elorder)+i];
        }
        elvec[nelem+nelemx+elbc] = new TElemento1d(upmatid,elorder,nodesbc);
    }
    
    
    for (int elbc=0; elbc<nelemy; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesyleft [elbc*(elorder)+i];
        }
        elvec[nelem+2*nelemx+elbc] = new TElemento1d(leftmatid,elorder,nodesbc);
        
    }
    
    
    
    for (int elbc=0; elbc<nelemy; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesyright [elbc*(elorder)+i];
        }
        elvec[nelem+2*nelemx+nelemy+elbc] = new TElemento1d(rightmatid,elorder,nodesbc);
        
    }
    
    int bctyped = 0;
    double contrstd = 0., contrhsd = 0.;
    mat = new TMaterialBC(downmatid,bctyped,contrstd,contrhsd);
    mesh.insertMaterial(mat);
    
    int bctypeu = 0;
    double contrstu = 0., contrhsu = 0.;
    mat = new TMaterialBC(upmatid,bctypeu,contrstu,contrhsu);
    mesh.insertMaterial(mat);
    
    int bctyper = 0;
    double contrstr = 0., contrhsr = 0.;
    mat = new TMaterialBC(rightmatid,bctyper,contrstr,contrhsr);
    mesh.insertMaterial(mat);
    
    int bctypel = 0;
    double contrstl = 0., contrhsl = 0.;
    mat = new TMaterialBC(leftmatid,bctypel,contrstl,contrhsl);
    mesh.insertMaterial(mat);
	
    //mesh.Print();
	
}


void CreateTestMeshTR(TMalha &mesh, int order)
{
    double xmin = 0.;
    double xmax = 1.;
    double ymin = 0.;
    double ymax = 1.;
    
    //int nelemx = (xmax-xmin)/h;
    //int nelemy = (ymax - ymin) / h;
    
	int nelemx = 1;
	int nelemy = nelemx;   //mismo numero de elementos
    
    int nelem = 2*nelemx*nelemy;
    
    
    
    int elorder = order;
    int nnodesx = nelemx*elorder + 1;  //nodes en x
    int nnodesy = nelemy*elorder + 1;  //mismo numero de nodes
	int nnodes = nnodesx*nnodesy;
    
	int matid = 1;
    double K = 1., B = 0., C[2] = { 0, 0 }, F(1.);  //elemento2d
    
	int leftmatid = -1;
    int rightmatid = -2;
	int downmatid = -3;
    int upmatid = -4;
    
    
    
    TVec<TNo> &nodevec = mesh.getNodeVec();
    nodevec.Resize(nnodes);
    
	TVecNum<double> coordx(nnodesx);
	TVecNum<double> coordy(nnodesy);
    
    for (int i=0; i<nnodesx; i++) {
        coordx[i] = xmin + (xmax-xmin)*i*1./(nnodesx-1);
	}
    
    for (int j=0; j<nnodesy; j++) {
        coordy[j] = ymin + (ymax-ymin)*j*1./(nnodesy - 1);
    }
    
    TVecNum<double> coord(2);
    
    for (int i = 0; i < nnodesx; i++)
    {
        for (int j = 0; j < nnodesy; j++)
        {
			coord[0] = coordx[j];
			coord[1] = coordy[i];
			nodevec[nnodesy*i+j].setData(coord);
        }
    }
    
    
    
    
    TVec<TElemento *> &elvec = mesh.getElementVec(); //inicializacao do objeto
    //elvec.Resize(nelem+4); //nelm elementos de interior + 2 elementos de BC;
    //elvec.Resize(nelem); //por lo pronto sin elementos de contorno;
	int nelem1D=nelemx*2+nelemy*2;
	elvec.Resize(nelem+nelem1D);
	
	
	TVec<int> nodes(1); //node index no interior do elemento;
	
    int count1d=0;
    
	for (int ely=0; ely<nelemy; ely++) {
		for (int elx = 0; elx<nelemx; elx++){
            for (int eli=0; eli<2; eli++) {
                switch (elorder){
                    case 1:
                        nodes.Resize(3);
                        
                        if (eli==0) {
                            
                            nodes[0]=0;					
                            nodes[1]=1;
                            nodes[2]=nnodesx;
                            
                            for (int i=0; i<nodes.Size(); i++)
                            {
                                nodes[i]= nodes[i]+elx+ely*nnodesx;
                                
                            }
                            
                        }else if(eli==1){
                            
                            nodes[0]=nnodesx+1;					
                            nodes[1]=nnodesx;
                            nodes[2]=1;
                            
                            for (int i=0; i<nodes.Size(); i++)
                            {
                                nodes[i]= nodes[i]+elx+ely*nnodesx;
                                
                            }
                            
                        }
                        
                        break;
                        
                    case 2:
                        nodes.Resize(6);
                        
                        if (eli==0) {
                            
                            nodes[0]=0;
                            nodes[1]=2;
                            nodes[2]=((nnodes-1)/nelemx)-2;
                            nodes[3]=1;
                            nodes[4]=nnodesx+1;
                            nodes[5]=nnodesx;
                            
                            for (int i=0; i<nodes.Size(); i++)
                            {
                                nodes[i]= nodes[i]+2*elx+2*ely*nnodesx;
                                
                            }
                            
                        }else if(eli==1){
                            
                            nodes[0]=(nnodes-1)/nelemx; 
                            nodes[1]=((nnodes-1)/nelemx)-2;
                            nodes[2]=2;
                            nodes[3]=nnodesx+nnodesy+1;
                            nodes[4]=nnodesx+1;
                            nodes[5]=nnodesx+2;
                            
                            for (int i=0; i<nodes.Size(); i++)
                            {
                                nodes[i]= nodes[i]+2*elx+2*ely*nnodesx;
                                
                            }
                            
                        }
                        
                        break;
                        
                }
               
                
                elvec[count1d++] = new TElementoTriangle(matid,elorder,nodes);
                
            }
            
        }
    }
	TMaterial *mat = new TMaterial2d(matid,K,C,B,F);
    mesh.insertMaterial(mat);
    
	
	/*TMatrix stiff, rhs;
     mesh.getElement(0)->CalcStiff(mesh,stiff,rhs);
     std::cout<<"matriz rigidez elemeto 0"<<std::endl;
     stiff.Print();
     */
	
	
	//condicoes de contorno
    
	TVecNum<double> nodesxdown(nnodesx);
    TVecNum<double> nodesxup(nnodesx);
	TVecNum<double> nodesyleft(nnodesy);
	TVecNum<double> nodesyright(nnodesy);
    
    for (int i=0; i<nnodesx; i++)
    {
        nodesxdown[i]= i;
        nodesyleft[i]=i*nnodesx;
        nodesxup[i]= i+nnodesx*(nnodesy-1);
        nodesyright[i]=i*nnodesx+(nnodesx-1);
    }
    
    
	
	
    TVec<int> nodesbc(elorder+1);
    
    for (int elbc=0; elbc<nelemx; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesxdown [elbc*(elorder)+i];
        }
        elvec[nelem+elbc] = new TElemento1d(downmatid,elorder,nodesbc);
    }
    
    for (int elbc=0; elbc<nelemx; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesxup [elbc*(elorder)+i];
        }
        elvec[nelem+nelemx+elbc] = new TElemento1d(upmatid,elorder,nodesbc);
    }
    
    
    for (int elbc=0; elbc<nelemy; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesyleft [elbc*(elorder)+i];
        }
        elvec[nelem+2*nelemx+elbc] = new TElemento1d(leftmatid,elorder,nodesbc);
        
    }
    
    
    
    for (int elbc=0; elbc<nelemy; elbc++) {
        
        for (int i=0; i<=elorder; i++) {
            nodesbc[i]=nodesyright [elbc*(elorder)+i];
        }
        elvec[nelem+2*nelemx+nelemy+elbc] = new TElemento1d(rightmatid,elorder,nodesbc);
        
    }
    
    int bctyped = 0;
    double contrstd = 0., contrhsd = 0.;
    mat = new TMaterialBC(downmatid,bctyped,contrstd,contrhsd);
    mesh.insertMaterial(mat);
    
    int bctypeu = 0;
    double contrstu = 0., contrhsu = 0.;
    mat = new TMaterialBC(upmatid,bctypeu,contrstu,contrhsu);
    mesh.insertMaterial(mat);
    
    int bctyper = 0;
    double contrstr = 0., contrhsr = 0.;
    mat = new TMaterialBC(rightmatid,bctyper,contrstr,contrhsr);
    mesh.insertMaterial(mat);
    
    int bctypel = 0;
    double contrstl = 0., contrhsl = 0.;
    mat = new TMaterialBC(leftmatid,bctypel,contrstl,contrhsl);
    mesh.insertMaterial(mat);
	
    //mesh.Print();
	
}


int TestMesh(TMalha *mesh)
{
    int res=0;
    cout << __PRETTY_FUNCTION__ <<endl;
    cout<<"Teste1: Verifica os valores da funcao de forma sob o elemento mestre"<<endl;
    int count=0;
    TVecNum<double> point(1,0), phi;
    TMatrix dphi;
    int testres=0;
    for (int order=1; order < 5; order++) {
        
        for (int pt=0; pt<order+1; pt++) {
            point[0]=pt*(2/(double)order)-1;
            
            mesh->getElement(0)->Shape1d(order,point, phi, dphi);
            
            if (phi[pt]==1) {
                testres++;
            }
        }
        
        if (testres == order+1) {
            cout <<"        ->Para ordem "<<order<<" os valores de phi est„o corretos."<<endl;
            count++;
        }
        testres=0;
        
    }
    
    if (count == 4) {
        cout <<"Teste1: Ok" <<endl;
        res++;
    }
    
    else cout <<"Teste1: Errado" <<endl;
    
    /**************************************************/
    count=0;
    point[0]=0;
    cout<<"Teste2: Verifica os valores da derivada da funcao de forma (dPhi) sob o elemento mestre."<<endl;
    
    for (int order=1; order<5; order++) {
        mesh->getElement(0)->Shape1d(order, point, phi, dphi);
        int testres=0;
        
        if (order==1) {
            if(fabs(dphi(0,0)+0.5)< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,1)-0.5)< 1.e-9)
            {
                testres++;
            }
            
            if (testres==dphi.Cols()) {
                cout <<"        ->Para ordem "<<order<<" os valores de dphi est„o corretos."<<endl;
                count++;
            }
        }
        
        if (order==2) {
            if(fabs(dphi(0,0)+0.5)< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,1))< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,2)-0.5)< 1.e-9)
            {
                testres++;
            }
            
            if (testres==dphi.Cols()) {
                cout <<"        ->Para ordem "<<order<<" os valores de dphi est„o corretos."<<endl;
                count++;
            }
            
        }
        
        if (order==3) {
            if(fabs(dphi(0,0)-0.0625)< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,1)+(1.6875))< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,2)-(1.6875))< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,3)+(0.0625))< 1.e-9)
            {
                testres++;
            }
            
            if (testres==dphi.Cols()) {
                cout <<"        ->Para ordem "<<order<<" os valores de dphi est„o corretos."<<endl;
                count++;
            }
            
        }
        
        if (order==4) {
            
            double val=0.16666666666667;
            if(fabs(dphi(0,0)-val)< 1.e-9)
            {
                testres++;
            }
            
            val = -1.33333333333333333333;
            if(fabs(dphi(0,1)-val)< 1.e-9)
            {
                testres++;
            }
            
            if(fabs(dphi(0,2))< 1.e-9)
            {
                testres++;
            }
            
            val = 1.33333333333333333333;
            if(fabs(dphi(0,3)-val)< 1.e-9)
            {
                testres++;
            }
            
            val = -0.16666666666667;
            if(fabs(dphi(0,4)-val)< 1.e-9)
            {
                testres++;
            }
            
            if (testres==dphi.Cols()) {
                cout <<"        ->Para ordem "<<order<<" os valores de dphi est„o corretos."<<endl;
                count++;
            }
            
        }
    }
    
    if (count == 4) {
        cout <<"Teste2: Ok" <<endl;
        res++;
    }
    
    else cout <<"Teste2: Errado" <<endl;
    
    /**************************************************/
    cout<<"Teste3: Verifica os valores do Jacobiano" <<endl;
    point.Zero();
    count=0;
    TMatrix jacobian, jacinv;
    double detjac;
    int NElement;
    NElement = mesh->getElementVec().Size()-2;
    
    TMatrix Cjacobian(1,1), Cjacinv(1,1);
    double Cdetjac;
    
    Cjacobian(0,0) = (1/(double)NElement)/2;
    Cjacinv(0,0) = 2/(1/(double)NElement);
    Cdetjac = Cjacobian(0,0);
    
    //Cjacobian.Print();
    //Cjacinv.Print();
    
    for(int iE=0; iE<NElement; iE++) {
        mesh->getElement(0)->Jacobian(point, jacobian, jacinv, detjac, *mesh);
        //(jacobian-Cjacobian).Print();
        //(jacinv-Cjacinv).Print();
        
        if (fabs(jacobian(0,0)-Cjacobian(0,0))<1.e-5 && fabs(jacinv(0,0)-Cjacinv(0,0))<1.e-5 && fabs(detjac-Cdetjac)<1.e-5) {
            count++;
        }
        
        
    }
    
    if (count == NElement) {
        cout <<"Teste3: Ok" <<endl;
        res++;
    }
    
    else cout <<"Teste3: Errado" <<endl;
    
    /**************************************************/
    cout<<"Teste4: Testa a matriz de rigidez para o elemento1 "<<endl;
    int order = mesh->getElement(0)->getNodeVec().Size()-1;
    
    count=0;
    point.Zero();
    phi.Zero();
    dphi.Zero();
    TMatrix dphix(dphi);
    
    TMatrix localStiff(order+1,order+1), localRhs(order+1,1,0);
    TMatrix ElementStiff1(order+order,order+order,0), ElementRhs1(order+order,1,0);
    
    TMaterial *mat = mesh->getMaterial(1);
    
    TIntRule1d  intrule(order+order);
    double weight;
    int IntPoints = intrule.NPoints();
    
    //Resultado dado pelo programa;
    for (int ip=0; ip<IntPoints; ip++) {
        
        intrule.Point(ip, point, weight);
        
        mesh->getElement(0)->Shape1d(order, point, phi, dphi);
        
        dphix = dphi;
        dphix.Zero();
        
        for(int i=0; i<order+1; i++)
        {
            dphix(0,i) = jacinv(0,0)*dphi(0,i);
        }
        
        weight *=fabs(detjac);
        
        mat->Contribute(weight, phi, dphix, localStiff, localRhs);
        std::cout<<localStiff(0,0)<<std::endl;
        
    }
    
    ElementRhs1 = localRhs;
    ElementStiff1 = localStiff;
    localStiff.Zero();
    localRhs.Zero();
    
    TMatrix ElementStiff2(order+order,order+order,0), ElementRhs2(order+order,1,0);
    
    //Resultado verdadeiro dado pelo test
    for (int ip=0; ip<IntPoints; ip++) {
        
        intrule.Point(ip, point, weight);
        
        mesh->getElement(0)->Shape(point, phi, dphi);
        
        for(int i=0; i<phi.Size(); i++)
        {
            dphix(0,i) = jacinv(0,0)*dphi(0,i);
        }
        
        weight *=fabs(detjac);
        
        double K = 1., B = -1., C=0., F(-1.);
        
        int i, j, nshape;
        nshape = phi.Size();
        
        for(i=0; i<nshape; i++)
        {
            localRhs(i,0) += weight*phi[i]*F;
            for(j=0; j<nshape; j++)
            {
                localStiff(i,j) += dphix(0,i)*dphix(0,j)*K*weight+
                phi[i]*dphix(0,j)*weight*C+
                phi[i]*phi[j]*B*weight;
            }
        }
        
        //std::cout << localStiff(0,0) << std::endl;
        
        
    }
    ElementRhs2 = localRhs;
    ElementStiff2 = localStiff;
    
    TMatrix Error(ElementStiff1);
    Error.Zero();
    
    Error = ElementStiff2 - ElementStiff1;
    // ElementStiff2.Print();
    //ElementStiff1.Print();
    //Error.Print();
    
    for (int i=0; i<Error.Cols(); i++) {
        for (int j=0; j<Error.Rows(); j++) {
            if (Error(i,j) < 0.001 || Error(i,j) > 0.001) {
                count++;
            }
        }
        
    }
    
    if (count == Error.Rows()*Error.Cols()) {
        res++;
        cout <<"Teste4: Ok\n";
    }
    else cout <<"Teste4: Errado\n";
    
    /**************************************************/
    cout<<"Teste5: Testa a assemblagem"<<endl;
    count=0;
    double BigNumber = 1.e12;
    
    //Resultado verdadeiro dado pelo teste
    
    int neq = mesh->getNodeVec().Size();
    
    TMatrix TestStiff(neq,neq,0.0);
    TMatrix TestRhs(neq,1);
    
    TMatrix ProgStiff(neq,neq,0.0);
    TMatrix ProgRhs(neq,1);
    
    for (int iE=0; iE<NElement; iE++) {
        
        TVec<int> nodes ;
        
        nodes = mesh->getElement(iE)->getNodeVec();
        
        int nnodes = nodes.Size();
        int in,jn;
        for(in=0; in<nnodes; in++)
        {
            TestRhs(nodes[in],0) += localRhs(in,0);
            for(jn=0; jn<nnodes; jn++)
            {
                TestStiff(nodes[in],nodes[jn]) += localStiff(in,jn);
            }
        }
    }
    
    //TestStiff.Print();
    
    TestStiff(0,0) += BigNumber;
    TestStiff(neq-1,neq-1) += BigNumber;
    
    TAnalysis Analysis(mesh);
    
    //Stiff pelo programa do usuario
    Analysis.Assemble(ProgStiff, ProgRhs);
    
    //ProgStiff.Print();
    
    
    
    TMatrix StiffError(TestStiff);
    StiffError.Zero();
    TMatrix RhsError(TestRhs);
    RhsError.Zero();
    
    StiffError = TestStiff-ProgStiff;
    //StiffError.Print();
    
    RhsError =TestRhs-ProgRhs;
    
    for (int i=0; i<StiffError.Rows(); i++) {
        
        if (fabs(RhsError(i,0))< 1.e-9){
            count++;
        }
        for (int j=0; j<StiffError.Cols(); j++) {
            if (fabs(StiffError(i,j))< 1.e-9) {
                count++;
            }
        }
    }
    
    if (count==neq*neq+neq){
        res++;
        cout<<"Test5: Ok\n";
    }
    
    else cout <<"Teste5: Errado\n";
    
    /**************************************************/
    cout<<"Teste6: Neuman BC"<<endl;
    count = 0;
    
    if (fabs(ProgStiff(0,0)-BigNumber-localStiff(0,0)) < 0.0001 || fabs(ProgStiff(neq-1,neq-1)-BigNumber-localStiff(0,order)) < 0.0001 ) {
        count++;
    }
    
    if (count == 1) {
        res++;
        cout <<"Teste6: Ok\n";
    }
    else cout <<"Teste6: Errado\n";
    
    cout << "---------------//----------------" <<endl;
    
    if (res==6) {
        res=1;
    }else{
        res=0;
    }
    
    return res;
}

//criacao de um ponteiro
int TestElement(TMalha *mesh)
{
    cout << __PRETTY_FUNCTION__ <<endl;
    
    int order = mesh->getElement(0)->getNodeVec().Size()-1;
    int res=0, count=0;
    int NElements = mesh->getElementVec().Size()-2;
    
    TVec<TNo> nodes = mesh->getNodeVec();
    
    cout<<"Teste1: Acesso a TVec<TElemento> em posicao inexistente."<<endl;
    int count1=0;
    
    for (int el=0; el < NElements; el++) {
        for (int index=0; index<=order; index++) {
            if (mesh->getElement(el)->getNodeVec().GetVal(index)==el*(order)+index){
                count1++;
            }
        }
    }
    
    if(count1==(NElements)*(order+1))
    {
        res=1;
    }
    if (res==1) {
        count++;
        cout<<"Teste1: Ok."<<endl;
    }
    if (res==0) {
        cout<<"Teste1: Errado."<<endl;
    }
    count1=0;   res=0;
    cout<<"Teste2: Testa a funcionalidade dos indices de TElemento com TNo."<<endl;
    for (int el=0; el<NElements; el++) {
        for (int index=0; index<=order; index++) {
            if (nodes[(order)*el+index].Co(0) == nodes[mesh->getElement(el)->getNodeVec().GetVal(index)].Co(0) ) {
                count1++;
            }
        }
    }
    if (count1 == (NElements)*(order+1)) {
        res=1;
    }
    if (res==1) {
        count++;
        cout<<"Teste2: Ok."<<endl;
    }
    if (res==0) {
        cout<<"Teste2: Errado."<<endl;
    }
    
    cout << "---------------//----------------" <<endl;
    res=0;
    if (count==2) {
        res=1;
    }
    return res;
}

int TestNodes(TMalha *mesh)
{
    cout << __PRETTY_FUNCTION__ <<endl;
    
    double xmin = 0.;
    double xmax = 1.;
    
    int order = mesh->getElement(0)->getNodeVec().Size()-1;
    int res=0, count=0;
    int NElements = mesh->getElementVec().Size()-2;
    int nnodes = NElements*order+1;
    
    TVec<TNo> nodes = mesh->getNodeVec();
    cout << "Teste1: Verificacao de acesso para node.Co[negativa e maior 1]." <<endl;
    try {
        nodes[0].Co(-1);
        nodes[1].Co(20);
    } catch (std::bad_exception) {
        res=1;
    }
    if (res==1) {
        count++;
        cout<<"Teste1: Ok."<<endl;
    }
    if (res==0) {
        cout<<"Teste1: Errado."<<endl;
    }
    
    res=0;
    cout <<"Teste2: Verificacao dos valores dos nos." <<endl;
    int count1=0;
    for (int i=0; i<nodes.Size(); i++) {
        double InstNode = xmin+(xmax-xmin)*i*1./(nnodes-1);
        
        if (nodes[i].Co(0)==InstNode) {
            count1++;
        }
    }
    if (count1==nodes.Size()) {
        res=1;
    }
    if (res==1) {
        count++;
        cout<<"Teste2: Ok."<<endl;
    }
    if (res==0) {
        cout<<"Teste2: Errado."<<endl;
    }
    
    res=0;
    cout <<"Teste3: Verificacao tamanho do vetor de nos." <<endl;
    if (nodes.Size()==nnodes) {
        res=1;
    }
    if (res==1) {
        count++;
        cout<<"Teste3: Ok."<<endl;
    }
    if (res==0) {
        cout<<"Teste3: Errado."<<endl;
    }
    
    cout << "---------------//----------------" <<endl;
    res=0;
    if (count==3) {
        res=1;
    }
    return res;
}


void TestOneDProblem(TMalha *mesh)
{
    int answer=0;
    
    answer += TestElement(mesh);
    answer += TestNodes(mesh);
    answer += TestMesh(mesh);
    
    if (answer==3) {
        cout <<"Seu teste passou no teste! :D"<<endl;
    }
    else{
        cout <<"Seu teste nao passou no teste! D:\nVerifique teste por teste e conserte os seus erros."<<endl;
    }
    
}
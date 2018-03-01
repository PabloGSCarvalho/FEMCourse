//
//  TestOneDProblem.cpp
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
        Os testes foram preparados com um proposito educacional,
        recomenda-se que o aluno entenda a funcionalidade de cada
        teste e posteriormente use com seu código caso a caso
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
#include "tmaterial1d.h"
#include "tmaterialbc.h"
#include "tanalysis.h"
#include "TIntRule1d.h"

using std::cout;
using std::endl;
using std::cin;

void CreateTestMesh(TMalha &mesh, int order);
int TestNodes(TMalha *mesh);
int TestElement(TMalha *mesh);
int TestMesh(TMalha *mesh);

void CreateTestMesh(TMalha &mesh, int order, double h);
void TestOneDProblem(TMalha *mesh);

void exact(TVec<double> &point,double &val, TVec<double> &deriv);

int main ()
{
    
    
    TMalha mesh;
    int order = 3;
    double h=1.0/8.0;
    CreateTestMesh(mesh,order,h);
    
    TAnalysis Analysis(&mesh);
    Analysis.Run();
    
    std::string filename("data.txt");
    Analysis.uh(filename);
    //exact(x, val, deriv);
    
    double energy=0.0,l2=0.0;
    
    Analysis.Error(exact, energy, l2);
    
    cout<<endl;
    cout<<energy<<endl;
    cout<<l2<<endl;
    

    //TestOneDProblem(&mesh);
    
    return 0;
}
void exact(TVec<double> &point,double &val, TVec<double> &deriv){


    double E=exp(1.0);
    TVec<double> x(1);
    x[0]=point[0];
    //x[0]=0;
    
//    val=(2.*(-1. + pow(E,sqrt(5.)*x[0]))*(sqrt(5.)*pow(E,10.*sqrt(5.)) + pow(E,20.*sqrt(5.)) - pow(E,sqrt(5.)*x[0]) + sqrt(5.)*pow(E,10.*sqrt(5.) + sqrt(5.)*x[0])))/
//    (pow(E,sqrt(5.)*x[0])*(1. + pow(E,20.*sqrt(5.))));
//    
//    deriv[0]=(2.*(5*pow(E,10.*sqrt(5.)) + sqrt(5.)*pow(E,20*sqrt(5.)) -
//                  sqrt(5.)*pow(E,2.*sqrt(5.)*x[0]) + 5.*pow(E,2*sqrt(5.)*(5. + x[0]))))/
//    (pow(E,sqrt(5.)*x[0])*(1. + pow(E,20.*sqrt(5.))));
    //cout<<val<<" ";
    //deriv.Print();
    val=(30. + 100.*pow(E,100.) - 130.*pow(E,10.*x[0]) - 3*x[0] + 3*pow(E,100.)*x[0])/(10.*(-1. + pow(E,100.)));
    deriv[0]=(-3. + 3*pow(E,100) - 1300*pow(E,10*x[0]))/(10.*(-1 + pow(E,100)));
    
    
}

void CreateTestMesh(TMalha &mesh, int order, double h)
{
    double xmin = 0.;
    double xmax = 10.;
    int nelem = (xmax-xmin)/h;
    int elorder = order;
    int nnodes = nelem*elorder+1;
    int matid = 1;
    double K = 1., B = 0., C=10., F(3.);
    int leftmatid = -1;
    int rightmatid = -2;
    
    TVec<TNo> &nodevec = mesh.getNodeVec();
    nodevec.Resize(nnodes);
    for (int i=0; i<nnodes; i++) {
        TVecNum<double> coord(2);
        coord[0] = xmin+(xmax-xmin)*i*1./(nnodes-1);
        nodevec[i].setData(coord);
    }
    
    TVec<TElemento *> &elvec = mesh.getElementVec(); //inicializacao do objeto
    elvec.Resize(nelem+2); //nelm elementos de interior + 2 elementos de BC;
    
    for (int el=0; el<nelem; el++) {
        TVec<int> nodes(elorder+1); //node index no interior do elemento;
        
        for (int i=0; i<=elorder; i++) {
            nodes[i]= el*(elorder)+i;
        }
        elvec[el] = new TElemento1d(matid,elorder,nodes);
    }
    
    {
        TVec<int> nodes(1);
        nodes[0] = 0;
        elvec[nelem] = new TElemento0d(leftmatid,0,nodes);
        nodes[0] = nnodes-1;
        elvec[nelem+1] = new TElemento0d(rightmatid,0,nodes);
    }

    
    TMaterial *mat = new TMaterial1d(matid,K,C,B,F);
    mesh.insertMaterial(mat);
    int bctype = 0;
    double contrst = 0., contrhs = 10.;
    mat = new TMaterialBC(leftmatid,bctype,contrst,contrhs);
    mesh.insertMaterial(mat);
    bctype = 0;
    contrst = 0., contrhs = 0.;
    mat = new TMaterialBC(rightmatid,bctype,contrst,contrhs);
    mesh.insertMaterial(mat);
}

int TestMesh(TMalha *mesh)
{
    int res=0;
    //cout << __PRETTY_FUNCTION__ <<endl;
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
            cout <<"        ->Para ordem "<<order<<" os valores de phi estão corretos."<<endl;
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
                cout <<"        ->Para ordem "<<order<<" os valores de dphi estão corretos."<<endl;
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
                cout <<"        ->Para ordem "<<order<<" os valores de dphi estão corretos."<<endl;
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
                cout <<"        ->Para ordem "<<order<<" os valores de dphi estão corretos."<<endl;
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
                cout <<"        ->Para ordem "<<order<<" os valores de dphi estão corretos."<<endl;
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
    double BigNumber = mesh->getMaterial(-1)->GetBigNumber();
    
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
    //cout << __PRETTY_FUNCTION__ <<endl;

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
        //cout << __PRETTY_FUNCTION__ <<endl;
        
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


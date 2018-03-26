

#include <iostream>
#include <math.h>
#include "TMatrix.h"
#include "TIntRule1d.h"
#include "TIntRuleQuad.h"
#include "TIntRuleTriangle.h"
#include "TVec.h"
#include "tpanic.h"
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

void test1D();

void testQuad();

void testTriangle();
//twodim Test
double Poli2D(TVecNum<double> &coord, int orderx, int ordery);
double Poli2DQuad(TVecNum<double> &coord, int orderx, int ordery);
void UXi(TVecNum<double> &coord, TVecNum<double> &uXi, TVecNum<double> &gradu);
TVecNum<double> X(TVecNum<double> &coordXi);
void Jacobian(TVecNum<double> &Coord, TMatrix &jacobian, double &detjac);

int main ()
{
//  test1D();
//  testTriangle();
//  testQuad();
    
    int order1 =19;
    TIntRuleQuad TestInt1(order1);
    
    //TestInt.Print(std::cout);
    

    double weight;
    TVecNum<double> CoordXi(2), uXi(1), gradu(2);
    
    double val1 = 0.;
    double val2x = 0.,val2y = 0.;
    int NPoint = TestInt1.NPoints();
    for (int i=0; i<NPoint; i++) {
        TestInt1.Point(i, CoordXi, weight);
        UXi(CoordXi,uXi,gradu);
        val1 = val1 + weight*uXi[0];
        val2x = val2x +weight*gradu[0];
        val2y = val2y +weight*gradu[1];
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1 << endl;
    cout<<"Resultado Integral Gradu(xi,eta):  "<< val2x << "," << val2y << endl;
    
    
    TMalha mesh;
    int Porder = 2;

    CoordXi.Zero();
    TMatrix jacobian, jacinv;
    double detjac;

    weight=0.;
    CoordXi.Zero(),uXi.Zero(),gradu.Zero();
    TIntRuleQuad TestInt2(order1);
    
    TVecNum<double> coord;
    double val3 = 0.,valA=0;
    double val4x = 0.,val4y = 0.;
    int NPoint2 = TestInt2.NPoints();
    TVec<int> nodes(NPoint2);
    
    TVec<TElemento *> &elvec = mesh.getElementVec();
    elvec.Resize(1);
    
    for (int i=0; i<NPoint2; i++) {
        nodes[i]=i;
    }
    elvec[0] = new TElementoQuad(1,1,nodes);
    
    for (int i=0; i<NPoint2; i++) {
        TestInt2.Point(i, CoordXi, weight);
        UXi(CoordXi,uXi,gradu);
        //coord = X(CoordXi);
        Jacobian(CoordXi, jacobian, detjac);
        valA = valA + detjac*weight;
        val3 = val3 + detjac*weight*uXi[0];
        val4x = val4x + detjac*weight*gradu[0];
        val4y = val4y + detjac*weight*gradu[1];
        
    }
    
    //
    
    
     cout<<"Resultado Área quad. mapeado:  "<< valA << endl;
     cout<<"Resultado Integral u(x,y):  "<< val3 << endl;
     cout<<"Resultado Integral Gradu(x,y):  "<< val4x << "," << val4y << endl;
    


    
  return 0;
  
}
//f[x_, y_] = x^2 + y^2 - 2 x y;
double Poli2D(TVecNum<double> &coord, int orderx, int ordery)
{
  return pow(coord[0],orderx) * pow(1.-coord[1],ordery);
}

double Poli2DQuad(TVecNum<double> &coord, int orderx, int ordery)
{
  return pow(coord[0],orderx) * pow(coord[1],ordery);
}

void UXi(TVecNum<double> &coord, TVecNum<double> &uxi, TVecNum<double> &gradu)
{
    uxi[0] = 3.+cos(3.*coord[1])*sin(4.*coord[0]);
    gradu[0] = 4.*cos(3.*coord[1])* cos(4.*coord[0]);
    gradu[1] = -3.*sin(3.*coord[1])* sin(4.*coord[0]);
}

TVecNum<double> X(TVecNum<double> &CoordXi)
{
    
    TVecNum<double> coord(2);

    coord[0] = 5.*CoordXi[0]+0.5*sin(3.*CoordXi[1]);
    coord[1] = 4.*CoordXi[1]+0.3*cos(10.*CoordXi[0]);

    return coord;
    
}

void Jacobian(TVecNum<double> &Coord, TMatrix &jacobian, double &detjac)
{

    jacobian.Resize(2, 2);
    jacobian.Zero();
    
    jacobian(0,0)=5.;
    jacobian(0,1)=1.5*cos(3.*Coord[1]);
    jacobian(1,0)=-3.*sin(10.*Coord[0]);
    jacobian(1,1)=4.;
    
    detjac = fabs(jacobian(0, 0)*jacobian(1, 1) - jacobian(1, 0)*jacobian(0, 1));
    
}


/*
 * Calcula os valores das funcoes de forma e suas derivadas
 * @point ponto onde calcular as funcoes de forma
 * @phi valores das funcoes de forma
 * @dphi valores das derivadas das funcoes de forma
 */

void Shape(TVec<double> &point, TVec<double> &phi, TMatrix &dphi, int order)
{
    
    if (order==1) {
        
        int Indices[2][2]={{0,3},{1,2}};
        
        TVec<double> coxi(1);
        coxi[0]=point[0];
        
        TVec<double> coeta(1);
        coeta[0]=point[1];
        
        TVec<double> phixi(2), phieta(2);
        TMatrix dphixi(1,2),dphieta(1,2);
        
        phi.Resize(4);
        dphi.Resize(2, 4);
        
        
        for (int xi=0; xi<order+1; xi++) {
            
            TElemento::Shape1d(order, coxi, phixi, dphixi);
            
            for (int eta=0; eta<order+1; eta++) {
                
                TElemento::Shape1d(order, coeta, phieta, dphieta);
                
                phi[Indices[xi][eta]]=phixi[xi]*phieta[eta];
                
                dphi(0,Indices[xi][eta])=dphixi(0,xi)*phieta[eta];
                dphi(1,Indices[xi][eta])=dphieta(0,eta)*phixi[xi];
                
            }
            
        }
        
    }
    
    if (order==2) {
        
        int Indices[3][3]={{0,7,3},{4,8,6},{1,5,2}};
        
        TVec<double> coxi(1);
        coxi[0]=point[0];
        
        TVec<double> coeta(1);
        coeta[0]=point[1];
        
        TVec<double> phixi(2), phieta(1);
        TMatrix dphixi(1,3),dphieta(1,3);
        
        phi.Resize(9);
        dphi.Resize(2, 9);
        
        
        for (int xi=0; xi<order+1; xi++) {
            
            TElemento::Shape1d(order, coxi, phixi, dphixi);
            
            for (int eta=0; eta<order+1; eta++) {
                
                TElemento::Shape1d(order, coeta, phieta, dphieta);
                
                phi[Indices[xi][eta]]=phixi[xi]*phieta[eta];
                
                dphi(0,Indices[xi][eta])=dphixi(0,xi)*phieta[eta];
                dphi(1,Indices[xi][eta])=dphieta(0,eta)*phixi[xi];
                
            }
            
        }
        
    }
    
}











double Poli1D(double coord, int polorder)
{
  double y;
  
  y = pow(coord, polorder);
  
  return y;
}

int TestTIntRule1d()
{
  cout <<__PRETTY_FUNCTION__<<endl;
  int result = 0;
  
  cout << "Teste1: Criacao de um TIntRule1d com ordem negativa" <<endl;
  try {
    TIntRule1d Integral(-1);
  } catch (std::bad_exception) {
    result = 1;
  }
  if (result == 0) {
    cout << "Teste1: Errado." <<endl;
  }
  else
    cout << "Teste1: Ok." <<endl;
  
  result = 0;
  cout << "Teste2: Criacao de um TIntRulde1d com ordem > 19" <<endl;
  try {
    TIntRule1d  Integral(20);
  } catch (std::bad_exception) {
    result = 1;
  }
  if (result == 0) {
    cout << "Teste2: Errado." <<endl;
  }
  else
    cout << "Teste2: Ok." <<endl;
  
  result = 0;
  cout << "Teste3: Verificacao se a ordem dada atribui corretamente o numero de pesos" <<endl;
  TIntRule1d  Integral(7);
  TIntRule1d  Integral2(18);
  
  if (Integral.NPoints() == 4) {
    result++;
  }
  else
    result = 0;
  
  if (Integral2.NPoints() == 10) {
    result++;
  }
  else
    result = 0;
  
  if (result == 2) {
    cout << "Teste3: Ok." <<endl;
  }
  else
    cout << "Teste3: Errado." <<endl;
  
  
  result = 0;
  cout << "Teste4: Verificacao dos pesos e coordenadas para uma ordem x." <<endl;
  TVecNum<double> testcoord(1);
  double peso;
  int pontoteste = 3;
  
  Integral.Point(pontoteste,testcoord,peso);
  
  if (testcoord[0] == 0.33998104358485626 && peso == 0.65214515486254609) {
    result++;
  }
  else
    result = 0;
  
  pontoteste = 8;
  Integral2.Point(pontoteste,testcoord,peso);
  
  if (testcoord[0] == -0.14887433898163122 && peso == 0.29552422471475287) {
    result++;
  }
  else
    result = 0;
  
  if (result != 2) {
    cout << "Teste1: Errado." <<endl;
  }
  if (result == 2) {
    cout << "Teste4: Ok." <<endl;
    result = 1;
  }
  
  cout << "---------------//----------------" <<endl;
  
  return result;
}

int TestPointFunctionFull()
{
  cout <<__PRETTY_FUNCTION__<<endl;
  int result = 0;
  
  cout << "Teste1: Teste completo de pesos e coordenadas para uma ordem x." <<endl;
  TIntRule1d  Teste(6);
  TVecNum<double> weights(4,0);
  TVec<double> valcoord(4.0);
  
  valcoord[0] = -0.86113631159405257;  weights[0] = 0.34785484513745385;
  valcoord[1] = 0.86113631159405257;   weights[1] = 0.34785484513745385;
  valcoord[2] = -0.33998104358485626;  weights[2] = 0.65214515486254609;
  valcoord[3] = 0.33998104358485626;   weights[3] = 0.65214515486254609;
  
  TVecNum<double> coord(1,0.);
  double peso;
  
  for(int i=0; i<Teste.NPoints() ;i++)
  {
    Teste.Point(i, coord, peso);
    
    if (coord[0] == valcoord[i] && peso == weights [i]) {
      result++;
    }
  }
  
  
  
  if (result != 4) {
    cout << "Teste1: Errado." <<endl;
  }
  if (result == 4) {
    cout << "Teste1: Ok." <<endl;
    result = 1;
  }
  
  cout << "---------------//----------------" <<endl;
  return result;
}

int TestPoliIntegration()
{
  cout <<__PRETTY_FUNCTION__<<endl;
  int result = 1;
  
  cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
  
  for (int order = 0; order <20; order++)
  {
    TIntRule1d TestInt(order);
    double val = 0;
    double weight;
    TVecNum<double> Coord(1);
    for (int polorder = 0; polorder <= order; polorder++)
    {
      val = 0.;
      for (int i=0; i<TestInt.NPoints(); i++) {
        TestInt.Point(i, Coord, weight);
        val = val + weight*Poli1D(Coord[0], polorder);
      }
      double exact = 1./(polorder+1.)*(1-pow(-1, polorder+1));
      double error = val-exact;
      
      if (fabs(error) >= 1.e-15) {
        std::cout << "A regra de integracao de ordem " << order << " nao integrou um polinome de ordem " << polorder << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
        result = 0;
      }
    }
    
  }
  if (result == 1) {
    cout << "Teste1: Ok." <<endl;
    result = 1;
  }
  else
    cout << "Teste1: Errado." <<endl;
  
  cout << "---------------//----------------" <<endl;
  return result;
}

int TestPoliIntegrationQuad()
{
  cout <<__PRETTY_FUNCTION__<<endl;
  int result = 1;
  
  cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
  
  for (int order = 0; order <20; order++)
  {
    TIntRuleQuad TestInt(order);
    
    //TestInt.Print(std::cout);
    
    double val = 0;
    double weight;
    TVecNum<double> Coord(2);
    for (int polorderx = 0; polorderx <= order; polorderx++)
    {
      for (int polordery=0; polordery <= order; polordery++)
      {
        val = 0.;
        int NPoint = TestInt.NPoints();
        for (int i=0; i<NPoint; i++) {
          TestInt.Point(i, Coord, weight);
          val = val + weight*Poli2DQuad(Coord, polorderx,polordery);
        }
        double exact = 1./(polorderx+1.)*(1-pow(-1, polorderx+1))*1./(polordery+1.)*(1-pow(-1, polordery+1));
        double error = val-exact;
        
        if (fabs(error) >= 1.e-14) {
          std::cout << "A regra de integracao quadrilatero de ordem " << order << " nao integrou um polinome de ordem " << polorderx << " , " << polordery << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
          result = 0;
        }
      }
    }
  }
  if (result == 1) {
    cout << "TesteQuad1: Ok." <<endl;
    result = 1;
  }
  else
    cout << "TesteQuad1: Errado." <<endl;
  
  cout << "---------------//----------------" <<endl;
  return result;
}

int TestPoliIntegrationTriangle()
{
  cout <<__PRETTY_FUNCTION__<<endl;
  int result = 1;
  
  cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
  
  for (int order = 0; order <20; order++)
  {
    TIntRuleTriangle TestInt(order);
    double val = 0;
    double weight;
    TVecNum<double> Coord(2);
    for (int polorderx = 0; polorderx <= order; polorderx++)
    {
      for (int polordery=0; polordery <= order-polorderx; polordery++)
      {
        val = 0.;
        int npoints = TestInt.NPoints();
        for (int i=0; i<npoints; i++) {
          TestInt.Point(i, Coord, weight);
          val = val + weight*Poli2D(Coord, polorderx,polordery);
        }
        double exact = 1./(1.+polorderx)/(2.+polorderx+polordery);
        double error = val-exact;
        
        if (fabs(error) >= 1.e-15) {
          std::cout << "A regra de integracao triangular de ordem " << order << " nao integrou um polinome de ordem " << polorderx << " , " << polordery << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
          result = 0;
        }
      }
    }
  }
  if (result == 1) {
    cout << "TesteTriangle1: Ok." <<endl;
    result = 1;
  }
  else
    cout << "TesteTriangle1: Errado." <<endl;
  
  cout << "---------------//----------------" <<endl;
  return result;
}


void test1D()
{
  int res = 0;
  int count = 0;
  
  res = TestTIntRule1d();
  if (res == 1) {
    count++;
  }
  
  res = TestPointFunctionFull();
  if (res == 1) {
    count++;
  }
  
  res = TestPoliIntegration();
  if (res == 1) {
    count++;
  }
  
  
  if (count == 1) {
    cout <<"\nSeu codigo passou no teste de Integracao para funcoes de uma dimensao. (:D)" << "\n---------------//----------------";
  }
  
}

void testQuad()
{
  int result = 0;
  result = TestPoliIntegrationQuad();
  if (result) {
    cout <<"\nSeu codigo passou no teste de Integracao sobre quadrilatero. (:D)" << "\n---------------//----------------";
    
  }
}

void testTriangle()
{
  int result = 0;
  result = TestPoliIntegrationTriangle();
  if (result) {
    cout <<"\nSeu codigo passou no teste de Integracao sobre triangulo. (:D)" << "\n---------------//----------------";
    
  }
}

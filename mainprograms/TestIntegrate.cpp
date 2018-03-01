

#include <iostream>
#include <math.h>
#include "TMatrix.h"
#include "TIntRule1d.h"
#include "TIntRuleQuad.h"
#include "TIntRuleTriangle.h"
#include "TVec.h"
#include "tpanic.h"

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

int main ()
{
  test1D();
  testTriangle();
  testQuad();
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
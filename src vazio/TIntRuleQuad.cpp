//
//  TIntRule1d.cpp
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//
//Integral cpp
#include "TIntRuleQuad.h"
#include "tpanic.h"
#include "TIntRule1d.h"


TIntRuleQuad::TIntRuleQuad(int order)
{
    SetOrder(order);
}

int TIntRuleQuad::NPoints()
{
    TIntRule1d Int1D(fOrder);
    
    return (Int1D.NPoints())*(Int1D.NPoints());
}

void TIntRuleQuad::Print(std::ostream &out)
{
    DebugStop();
}

void TIntRuleQuad::SetOrder(int order)
{
    if (order<0||order>19) {
        DebugStop();
    }
    
    fOrder=order;
    
}

void TIntRuleQuad::Point(int p, TVec<double> &co, double &weight)
{
    if(p<0||p>=NPoints()){
        DebugStop();
    }
    
    TIntRule1d Int1Dx(fOrder);
    TIntRule1d Int1Dy(fOrder);
    
    fPoints.Resize(NPoints(), 2);
    fWeights.Resize(NPoints());

    for (int i=0; i<Int1Dx.NPoints(); i++) {
        
        Int1Dx.Point(i, co, weight);
        TVecNum<double> coX(1);
        double weightX;
        coX[0]=co[0];
        weightX=weight;
        
        for (int j=0; j<Int1Dy.NPoints(); j++) {
            
            Int1Dy.Point(j, co, weight);
        
            fPoints(j+i*Int1Dy.NPoints(),0)=co[0];
            fPoints(j+i*Int1Dy.NPoints(),1)=coX[0];
            
            fWeights[j+i*Int1Dy.NPoints()]=weightX*weight;
        }
        
    }
    
    co.Resize(2);
    
    co[0]=fPoints(p,0);
    co[1]=fPoints(p,1);
    
    weight=fWeights[p];
    
    
}

void TIntRuleQuad::gaulegQuad(const double x1, const double x2, TVecNum<double> &x, TVecNum<double> &w)
{
    
    TIntRule1d IntGauss1Dx(fOrder);
    TIntRule1d IntGauss1Dy(fOrder);
    double nPoints = x.Size();
    TVecNum<double> weightx(x.Size()), coX(nPoints);
    TVecNum<double> weighty(x.Size()), coY(nPoints);
    
    IntGauss1Dx.gauleg(x1, x2, coX, weightx);
    IntGauss1Dy.gauleg(x1, x2, coY, weighty);
    
    x.Resize(2*nPoints*nPoints);
    w.Resize(nPoints*nPoints);
    
    for (int i = 0; i<nPoints; i++) {
        
        for (int j = 0; j<nPoints; j++) {
            w[j+i*nPoints]=weightx[j]*weighty[i];
            x[j+i*nPoints]=coX[j];
            x[j+i*nPoints+nPoints*nPoints]=coY[i];
        }
    }
    
}





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

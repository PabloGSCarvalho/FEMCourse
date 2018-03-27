//
//  TIntRule1d.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__TIntRuleQuad__
#define __FemSC__TIntRuleQuad__

#include <stdio.h>
#include "TVec.h"
#include "TMatrix.h"

class TIntRuleQuad
{

  int fOrder;
  
  TMatrix fPoints;
  
  TVec<double> fWeights;
  

public:
  
  TIntRuleQuad();
  
    TIntRuleQuad(int order);
  
  void SetOrder(int order);
    
    int NPoints();
    
    void Point(int p, TVec <double> &co, double &weight);
    
    void gaulegQuad(const double x1, const double x2, TVecNum<double> &x, TVecNum<double> &w);
    
    void Print(std::ostream &out);
};


#endif /* defined(__FemSC__TIntRule1d__) */

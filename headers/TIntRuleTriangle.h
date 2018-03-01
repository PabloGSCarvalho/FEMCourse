//
//  TIntRule1d.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__TIntRuleTriangle__
#define __FemSC__TIntRuleTriangle__

#include <stdio.h>
#include "TVec.h"
#include "TMatrix.h"

class TIntRuleTriangle
{
  int fOrder;
  
  TMatrix fPoints;
  
  TVec<double> fWeights;
  
  
public:
  
  TIntRuleTriangle();
  
  TIntRuleTriangle(int order);
  
  void SetOrder(int order);
  
  int NPoints();
  
  void Point(int p, TVec <double> &co, double &weight);
  
};


#endif /* defined(__FemSC__TIntRule1d__) */

//
//  TIntRule1d.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__TIntRule1d__
#define __FemSC__TIntRule1d__

#include <cmath>
#include <stdio.h>
#include "TVecNum.h"

class TIntRule1d
{
  
  int fOrder;
  
  TVecNum<double> fPoints;
  
  TVecNum<double> fWeights;
    
public:
  
    TIntRule1d();
    
    TIntRule1d(int order);
    
    void SetOrder(int order);
    
    int NPoints();
    
    void Point(int p, TVec<double> &co, double &weight);
    
    void gauleg(const double x1, const double x2, TVecNum<double> &x, TVecNum<double> &w);
    
};


#endif /* defined(__FemSC__TIntRule1d__) */

//
//  IntRule.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRule__
#define __FemSC__IntRule__

#include <cmath>
#include <stdio.h>
#include "TVecNum.h"

class IntRule
{
  
protected:
    
    int fOrder;

    TMatrix fPoints;
    
    TVec<double> fWeights;

public:
  
    //Default Constructor of integration rule
    IntRule();
    
    
    IntRule(int order);
    
    ~IntRule();
    
    virtual void operator=(const IntRule &copy);
    
    IntRule(const IntRule &copy);
    
    virtual void SetOrder(int order)
    {
        fOrder=order;
        
    }
    
    virtual int NPoints() const;
    
    virtual void Point(int p, TVec<double> &co, double &weight) const;
    
    void Print(std::ostream &out) const;
    
};


#endif /* defined(__FemSC__TIntRule1d__) */

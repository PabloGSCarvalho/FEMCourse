//
//  IntRuleTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRuleTriangle__
#define __FemSC__IntRuleTriangle__

#include <stdio.h>
#include "TVec.h"
#include "TMatrix.h"

class IntRuleTriangle : public IntRule
{
  
public:
  
  IntRuleTriangle();
  
  IntRuleTriangle(int order);
  
  virtual void SetOrder(int order);
  
};


#endif /* defined(__FemSC__TIntRule1d__) */

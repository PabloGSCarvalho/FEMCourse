//
//  IntRuleTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRuleTetrahedron__
#define __FemSC__IntRuleTetrahedron__

#include <stdio.h>
#include "TVec.h"
#include "TMatrix.h"

class IntRuleTetrahedron : public IntRule
{
  
    
public:
  
  IntRuleTetrahedron();
  
  IntRuleTetrahedron(int order);
  
  virtual void SetOrder(int order);
    
};


#endif /* defined(__FemSC__TIntRule1d__) */

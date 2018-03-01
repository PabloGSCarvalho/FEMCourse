//
//  Header.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef FemSC_Panic_h
#define FemSC_Panic_h

#include <iostream>
#include <exception>

static void DebugStop()
{
    std::cout << "Your chance to put a breakpoint here\n";
    std::bad_exception myex;
    throw myex;

}
#endif

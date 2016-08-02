#pragma once
#ifndef RSSERROR_H
#define RSSERROR_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <map>
#include "GenomeCopyNumber.h"
#include "ChrCopyNumber.h"


class RSSerror
{
    public:
        RSSerror();
};

long double calculateRSS(GenomeCopyNumber & samplecopynumber, int ploidy);

#endif // RSSERROR_H

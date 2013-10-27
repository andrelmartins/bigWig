//
//  bwgExtra.h
//  bigWigXC
//
//  Created by Andre Martins on 10/27/13.
//  Copyright (c) 2013 Andre Martins. All rights reserved.
//

#ifndef bigWigXC_bwgExtra_h
#define bigWigXC_bwgExtra_h

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sig.h"
#include "obscure.h"
#include "dystring.h"
#include "bPlusTree.h"
#include "cirTree.h"
#include "rangeTree.h"
#include "udc.h"
#include "zlibFace.h"
#include "bbiFile.h"
#include "bwgInternal.h"
#include "bigWig.h"

struct bbiInterval *bigWigIntervalQueryNoClip(struct bbiFile *bwf, char *chrom, bits32
                                              start, bits32 end,
                                              struct lm *lm);

#endif

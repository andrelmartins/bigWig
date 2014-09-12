#include "bwgExtra.h"

/* Adapted from Kent source (included jkweb folder) to not perform clipping */

struct bbiInterval *bigWigIntervalQueryNoClip(struct bbiFile *bwf, char *chrom, bits32
                                        start, bits32 end,
                                        struct lm *lm)
/* Get data for interval.  Return list allocated out of lm. */
{
  if (bwf->typeSig != bigWigSig)
    errAbort("Trying to do bigWigIntervalQuery on a non big-wig file.");
  bbiAttachUnzoomedCir(bwf);
  struct bbiInterval *el, *list = NULL;
  struct fileOffsetSize *blockList = bbiOverlappingBlocks(bwf, bwf->unzoomedCir,
                                                          chrom, start, end, NULL);
  struct fileOffsetSize *block, *beforeGap, *afterGap;
  struct udcFile *udc = bwf->udc;
  boolean isSwapped = bwf->isSwapped;
  float val;
  int i;
  
  
  /* Set up for uncompression optionally. */
  char *uncompressBuf = NULL;
  if (bwf->uncompressBufSize > 0)
    uncompressBuf = needLargeMem(bwf->uncompressBufSize);
  
  /* This loop is a little complicated because we merge the read requests for efficiency, but we
   * have to then go back through the data one unmerged block at a time. */
  for (block = blockList; block != NULL; )
  {
    /* Find contigious blocks and read them into mergedBuf. */
    fileOffsetSizeFindGap(block, &beforeGap, &afterGap);
    bits64 mergedOffset = block->offset;
    bits64 mergedSize = beforeGap->offset + beforeGap->size - mergedOffset;
    udcSeek(udc, mergedOffset);
    char *mergedBuf = needLargeMem(mergedSize);
    udcMustRead(udc, mergedBuf, mergedSize);
    char *blockBuf = mergedBuf;
    
    /* Loop through individual blocks within merged section. */
    for (;block != afterGap; block = block->next)
    {
      /* Uncompress if necessary. */
      char *blockPt, *blockEnd;
      if (uncompressBuf)
      {
        blockPt = uncompressBuf;
        int uncSize = zUncompress(blockBuf, block->size, uncompressBuf, bwf->uncompressBufSize);
        blockEnd = blockPt + uncSize;
      }
      else
      {
        blockPt = blockBuf;
        blockEnd = blockPt + block->size;
      }
      
      /* Deal with insides of block. */
      struct bwgSectionHead head;
      bwgSectionHeadFromMem(&blockPt, &head, isSwapped);
      switch (head.type)
      {
        case bwgTypeBedGraph:
        {
          for (i=0; i<head.itemCount; ++i)
          {
            bits32 s = memReadBits32(&blockPt, isSwapped);
            bits32 e = memReadBits32(&blockPt, isSwapped);
            val = memReadFloat(&blockPt, isSwapped);
            if (s < end && e > start)
            {
              lmAllocVar(lm, el);
              el->start = s;
              el->end = e;
              el->val = val;
              slAddHead(&list, el);
            }
          }
          break;
        }
        case bwgTypeVariableStep:
        {
          for (i=0; i<head.itemCount; ++i)
          {
            bits32 s = memReadBits32(&blockPt, isSwapped);
            bits32 e = s + head.itemSpan;
            val = memReadFloat(&blockPt, isSwapped);
            if (s < end && e > start)
            {
              lmAllocVar(lm, el);
              el->start = s;
              el->end = e;
              el->val = val;
              slAddHead(&list, el);
            }
          }
          break;
        }
        case bwgTypeFixedStep:
        {
          bits32 s = head.start;
          bits32 e = s + head.itemSpan;
          for (i=0; i<head.itemCount; ++i)
          {
            val = memReadFloat(&blockPt, isSwapped);
            bits32 clippedS = s, clippedE = e;
            if (clippedS < end && clippedE > start)
            {
              lmAllocVar(lm, el);
              el->start = clippedS;
              el->end = clippedE;
              el->val = val;
              slAddHead(&list, el);
            }
            s += head.itemStep;
            e += head.itemStep;
          }
          break;
        }
        default:
          internalErr();
          break;
      }
      assert(blockPt == blockEnd);
      blockBuf += block->size;
    }
    freeMem(mergedBuf);
  }
  freeMem(uncompressBuf);
  slFreeList(&blockList);
  slReverse(&list);
  return list;
}

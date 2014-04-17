ifeq ($(KENTHOME),)
	KENTHOME :=/home/cgd24/src/jksrc/kent
endif

include $(KENTHOME)/src/inc/common.mk

ifeq (${MACHTYPE},)
	MACHTYPE=x86_64
endif

L += ${SOCKETLIB} -lm -lz
#MYLIBDIR = ./lib/${MACHTYPE}
#MYLIBS = ${MYLIBDIR}/jkweb.a
MYLIBS = jkweb.a

# build our MAKEFLAGS
# R needs spaces to be escaped ...
#
MFLGS := "CFLAGS=-I./inc ${COPT} ${CFLAGS}"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

#
# rebuild jkweb with -fPIC
#
O = aliType.o annoFilter.o annoFormatter.o annoFormatTab.o \
    annoGrator.o annoGratorQuery.o annoRow.o annoStreamer.o annoStreamVcf.o \
    apacheLog.o asParse.o axt.o axtAffine.o bamFile.o base64.o \
    basicBed.o bbiRead.o bbiWrite.o bigBed.o binRange.o bits.o \
    blastOut.o blastParse.o boxClump.o boxLump.o bPlusTree.o \
    bwgCreate.o bwgQuery.o bwgValsOnChrom.o \
    cda.o chain.o chainBlock.o chainConnect.o chainToAxt.o chainToPsl.o \
    cheapcgi.o cirTree.o codebias.o colHash.o common.o correlate.o crTree.o \
    dgRange.o diGraph.o dlist.o dnaLoad.o dnaMarkov.o dnaMotif.o dnaseq.o \
    dnautil.o dtdParse.o dystring.o \
    emblParse.o errCatch.o errabort.o \
    fa.o ffAli.o ffScore.o filePath.o fixColor.o flydna.o fof.o \
    font/mgCourier10.o font/mgCourier12.o font/mgCourier14.o font/mgCourier18.o \
    font/mgCourier24.o font/mgCourier34.o font/mgCourier8.o font/mgHelvetica10.o \
    font/mgHelvetica12.o font/mgHelvetica14.o font/mgHelvetica18.o font/mgHelvetica24.o \
    font/mgHelvetica34.o font/mgHelvetica8.o font/mgHelveticaBold10.o font/mgHelveticaBold12.o \
    font/mgHelveticaBold14.o font/mgHelveticaBold18.o font/mgHelveticaBold24.o \
    font/mgHelveticaBold34.o font/mgHelveticaBold8.o font/mgSixhi6.o font/mgSail8.o \
    font/mgTimes10.o font/mgTimes12.o font/mgTimes14.o font/mgTimes18.o \
    font/mgTimes24.o font/mgTimes34.o font/mgTimes8.o font/mgMenlo12.o \
    fuzzyShow.o \
    gapCalc.o gdf.o gemfont.o genomeRangeTree.o \
    gfNet.o gff.o gff3.o gfxPoly.o gifLabel.o \
    hacTree.o hash.o histogram.o hmmPfamParse.o hmmstats.o htmlPage.o htmshell.o \
    https.o intExp.o intValTree.o internet.o itsa.o \
    jointalign.o jpegSize.o keys.o knetUdc.o kxTok.o linefile.o lineFileOnBigBed.o localmem.o log.o \
    maf.o mafFromAxt.o mafScore.o md5.o memalloc.o memgfx.o metaWig.o mgCircle.o \
    mgPolygon.o mime.o net.o nib.o nibTwo.o nt4.o numObscure.o \
    obscure.o oldGff.o oligoTm.o options.o osunix.o pairHmm.o peakCluster.o \
    phyloTree.o pipeline.o portimpl.o pngwrite.o psGfx.o psPoly.o pscmGfx.o \
    psl.o pslGenoShow.o pslShow.o pslTbl.o pslTransMap.o pthreadWrap.o \
    qa.o quickHeap.o quotedP.o \
    ra.o rainbow.o rbTree.o rangeTree.o regexHelper.o repMask.o \
    rle.o rnautil.o rqlEval.o rqlParse.o rudp.o \
    scoreWindow.o seg.o seqOut.o seqStats.o servBrcMcw.o servCrunx.o \
    servcis.o servcl.o servmsII.o servpws.o shaRes.o slog.o snof.o \
    snofmake.o snofsig.o spaceSaver.o spacedColumn.o spacedSeed.o \
    splatAli.o sqlList.o sqlNum.o subText.o sufa.o sufx.o synQueue.o \
    tabRow.o textOut.o tokenizer.o trix.o twoBit.o \
    udc.o vcf.o vGfx.o vPng.o verbose.o \
    wildcmp.o wormdna.o \
    xAli.o xa.o xap.o xenshow.o xmlEscape.o xp.o zlibFace.o

JKWEBO := $(patsubst %.o,jkweb/%.o,$(O))

all: prepare bigwig.so

prepare:
	ln -f -s $(KENTHOME)/src/inc
	ln -f -s $(KENTHOME)/src/lib
	[ -d jkweb ] || mkdir -p jkweb/font

jkweb.a:
	@for obj in $(O); do echo "CC $${obj}"; \
	${CC} ${COPT} -I./inc ${CFLAGS} ${HG_DEFS} ${LOWELAB_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -fPIC -o jkweb/$${obj} -c lib/`dirname $${obj}`/`basename $${obj} o`c; done
	ar rcus jkweb.a $(JKWEBO)

#jkweb/%.o: lib/%.c
#	${CC} ${COPT} -I./inc ${CFLAGS} ${HG_DEFS} ${LOWELAB_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -fPIC -o $@ -c $<

bigwig.so: prepare jkweb.a bigWig_R.c utils_R.c bwgExtra.h bwgExtra.c bw_query.h bw_query.c bw_base.h bw_base.c bigWig_R2.c
	MAKEFLAGS=${MFLGSESC} R CMD SHLIB -o $@ bigWig_R.c utils_R.c bwgExtra.c  bw_base.c bw_query.c bigWig_R2.c ${MYLIBS} ${L}
	rm inc lib

clean:
	rm -f bigwig.so *.o lib inc jkweb/*.o jkweb/font/*.o jkweb.a

.PHONY: prepare

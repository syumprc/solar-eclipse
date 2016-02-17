all:	mlink

# Makefile for FASTLINK 3.0P
# Chris Hyams (cgh@cs.rice.edu)
#
# History:
# 5/23/95 cgh -- created
# 6/7/95 cgh -- added comments for v2.3P beta release
# 7/17/95 cgh -- updated for portability to most versions of make
#                consolidated Makefile.common
# 7/18/95 cgh -- added VMS-specific section written by
#                Kimmo Kallio (Kimmo.Kallio@ktl.fi)
# 12/95   aas -- eliminated MAXWORKINGSET and maxhap
# 2/97    Anna Becker/aas -- added llopbrk.c for use in unknown
#

# Please read the file README.makefile that came
# with your FASTLINK distribution for an overview
# on using this Makefile

# If you plan to run parallel FASTLINK, please refer
# to README.parallel for more details.

# One note:
# Most versions of make allows you to override Makefile variables
# from the command line.  This Makefile includes a variable EXTRAFLAGS
# that is intentionally left empty to allow you to add command-line
# flags without editing the Makefile.  EXTRAFLAGS is automatically
# passed to each of the different targets.
#
# For example, if you wanted to pass along an additional include
# directive to the compiler when making ILINK, you could type:
#
#  make ilink EXTRAFLAGS="-I/usr/ucb/include -I/usr/lib/include"
#
# Note that if you want to add more than one flag, you need to protect
# the argument in double quotes (") as above.


######################################################
################## Miscellaneous #####################
######################################################
# The following variables apply to compilation of
# *both* sequential *and* parallel FASTLINK

# Destination bin directory for executables.  Note: this
# directory *must* exist, or compilation will fail during
# linkage.
BINDIR = ./

# Compiler optimization flag.  Comment this out
# if your compiler 
OPTFLAG = -O

# If you use the cc compiler, uncomment this flag
# to use K&R style function prototypes.
# Note: FASTLINK will *not* compile under cc without
# this flag set.
#PROTO = -DKNR_PROTO

# System-dependent compilation flags.
# - Uncomment -DSOLARIS here if you are compiling on a
#   Sun Solaris platform.
# - If you are using gcc under SunOS 4.1.x, you can
#   use -DSUNOS_4_1_X_GCC to get rid of compiler
#   warnings.
# - Uncomment -D_POSIX_SOURCE if you are compiling for AIX.
# - If you have other system-specific compiler flags that you
#   need to pass to FASTLINK, uncomment the line "SYSDEP =",
#   and add the flag there.
# Otherwise, leave this alone.
#SYSDEP = -DSUNOS_4_1_X_GCC
#SYSDEP = -DSOLARIS
#SYSDEP = -D_POSIX_SOURCE
#SYSDEP =


######################################################
############### Sequential Code Variables ############
######################################################
# Edit the following when compiling sequential FASTLINK.

# C compiler for sequential code
CC = gcc

# Uncomment this if you are compiling for DOS, or if
# you want to turn off checkpointing in FASTLINK.
#CKPT = -DDOS


######################################################
############ Parallel Code Variables #################
######################################################
# These variables apply only to parallel FASTLINK.

# C compiler for parallel code
PARCC = cc

# Parallel support library.
# If you are using TreadMarks, leave these
# commented out.
# If you are using P4 or PARMACS, uncomment the
# appropriate variable.  Make sure that no more than
# one of the following lines is uncommented!
#PARLIB = -DIS_PARMACS=1
PARLIB = -DIS_P4=1

# Parallel include path.
# If you are using TreadMarks or P4, uncomment the
# appropriate variable.  If you are using PARMACS,
# comment these out, or supply your own, appropriate path.
# Make sure that no more than one of these variables
# is uncommented!
#PARINCLPATH = -I$(TmkDIR)/include
#PARINCLPATH = -I$(P4_INCLDIR)


########### TreadMarks specifics ###############

# Set this to the root directory of your TreadMarks installation
TmkDIR = /net/blitzen/Tmk-0.9.6.FASTLINK

# Target architecture.
# For example, if your TreadMarks library is in
# $(TmkDIR)/lib.mips, you should set this variable
# to mips.  Consult your TreadMarks documentation
# for more details.
ARCH = sparc

# If you are using TreadMarks on an ATM network,
# set this to the location of libatm.a on your
# system.
ATMDIR = /usr/fore/lib


########### p4 specifics #######################

# Set this to the root of your p4 installation
P4_HOME_DIR = /usr/local/lib/p4

# Uncomment this if you are running in parallel on SOLARIS
# This is *very* important -- FASTLINK will not compile
# without it.
#P4_MDEP_LD = -lsocket -lnsl -lthread


######### Parallel C Compiler Flags #############
# Edit the following variables when compiling
# parallel FASTLINK to suit your particular run.
#
# There are 3 sets of variables provided -- one each for
# ILINK, MLINK, and LINKAP.  The variables are 
# **PARMEM, **PRECOMP, where ** represents one of IL, ML,
# or LI, for ILINK, MLINK and LINKMAP respectively.
# They are briefly described below.  Please refer to the
# files README.constants, README.memory, README.parallel,
# and README.makefile for more details.
#
# **PARMEM
#    Defines the variable LESSMEMORY, instructing FASTLINK
#    to use an alternate algorithm which requires less memory,
#    but is slower.  Please refer to the README.memory
#    and README.makefile for details.
#    Uncomment this if you wish to use less memory.
#    Note: setting this value to 0 *or* 1 will both have the
#    same effect of defining LESSMEMORY.  If you don't want
#    less memory, be sure to comment this out.
# **PRECOMPUTE
#    Defines the variable PRECOMPUTE as 1 or 0.  PRECOMPUTE=1
#    instructs FASTLINK to use an alternate algorithm which
#    requires more memory, and is faster.  PRECOMPUTE=0
#    uses the normal approach.  Please refer to the README.makefile,
#    and README.parallel for details.
#    Unlike **PARMEM, this value must be set to 0 or 1.

### ILINK flags -- used when ILFLAGS is chosen
# edit these as you wish
#ILPARMEM	= -DLESSMEMORY
ILPRECOMP	= -DPRECOMPUTE=1

### MLINK flags -- used when MLFLAGS is chosen
# edit these as you wish
#MLPARMEM	= -DLESSMEMORY
MLPRECOMP	= -DPRECOMPUTE=1


## LINKMAP flags -- used when LIFLAGS is chosen
# edit these as you wish

#LIPARMEM	= -DLESSMEMORY
LIPRECOMP	= -DPRECOMPUTE=1


######################################################
############ Unix Utilities on Your System ###########
######################################################
# If any of the following utilities are different on
# your system, you can change them here.

# RM is used for "make clean", so if you want to
# use that, you'll need something defined here.
RM = rm -f


######################################################
################## VMS Specific ######################
######################################################
# 18-Jul-1995 Kimmo Kallio (Kimmo.Kallio@ktl.fi)
# Uncomment following lines if compiling on VMS 
# system, up to (but not including) the line that starts
# "END: VMS".
#
#CC     = mcr sys$disk:[]maketool_vms.
#BINDIR = 
#%ifndef CKPT
##
## If CKPT is not defined, use empty value
##
#CKPT = 
#%endif
#
## Unix C-compiler environment, not needed on VMS
#PROTO = 
#SYSDEP =
#CFLAGS =
#SMALLCFLAGS =
#OPTFLAG =
#LDFLAGS =
#DEBUGFLAG = 
#WARNFLAG =
#BAR_OUT =
#LOAD_OUT =
#
## Parallel stuff, not in use on VMS
#PARCC = 
#PARLIB = 
#PARINCLPATH = 
#P4_INCLDIR =
#P4_MDEP_LD =
#MDEP_LIBS =
#PARDEBUG =
#PARTHETAS =
#ILPARMEM =
#MLPARMEM =
#LIPARMEM =
#DO_LOAD =
#
#vmsfast: maketool_vms installfast 
#
#vmsslow: maketool_vms installslow 
#
#vmscleanexe: 
#   delete/log ilink.;*
#   delete/log mlink.;*
#   delete/log linkmap.;*
#   delete/log lodscore.;*
#
#vmsclean:
#   purge/log *.o,*./exclude=(makefile.)
#
#maketool_vms: maketool_vms.c
#  @maketool_vms "-O maketool_vms maketool_vms.c"
#
# END: VMS -- uncomment lines above here for VMS
#


############################################################
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
############################################################
#                                                          #
# You shouldn't need to change anything beyond this point. #
#                                                          #
############################################################
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
############################################################


######################################################
################## miscellaneous #####################
######################################################

# "Slow" compiler flag
SFLAG = -DLESSMEMORY

# Extra flags passed to each target.  Do not change this!
EXTRAFLAGS = 


######################################################
############### sequential code variables ############
######################################################

# For debugging
#DEBUGFLAG = -g

# Compiler flags for sequential code
CFLAGS = \
	 $(CKPT) \
	$(DEBUGFLAG) \
	$(EXTRAFLAGS) \
	$(OPTFLAG) \
	$(PROTO) \
	$(SYSDEP) \
	$(WARNFLAG)

# For sequential programs that don't take -DDOS, etc.
SMALLCFLAGS = \
	 $(DEBUGFLAG) \
	$(EXTRAFLAGS) \
	$(OPTFLAG) \
	$(PROTO) \
	$(SYSDEP) \
	$(WARNFLAG)

# Linker flags
LDFLAGS =  -lm

######################################################
############ parallel code variables #################
######################################################

################# TreadMarks stuff ###################

# TreadMarks library location
TmkLIB	= $(TmkDIR)/lib.$(ARCH)

# TreadMarks profiling directories
PROFDIR = /home/rrk/Ymk
PROFLIB = $(PROFDIR)/$(ARCH)
LIBPERF = /home/rrk/Ymk/gdb-4.13


################# P4 stuff ###########################

# p4 include path, library directory, and libraries
P4_INCLDIR = $(P4_HOME_DIR)/include
P4_LIBDIR = $(P4_HOME_DIR)/lib
P4_LIBS = $(P4_LIBDIR)/libp4.a $(MDEP_LIBS)

################# linker flags #######################

# LDFLAGS for each library
ATMLD   = $(LDFLAGS) -L$(TmkLIB) -lTmk_atm -L$(ATMDIR) -latm
UDPLD   = $(LDFLAGS) -L$(TmkLIB) -lTmk_udp
P4LD    = $(LDFLAGS) -L$(P4_LIBDIR) -lp4 $(P4_MDEP_LD)

########### Various compiler flags ###################

# Parallel flag
DOPAR = -DPARALLEL=1

# Supress checkpointing
NOCKPT = -DDOS

# Supress Barrier output
#BAR_OUT = -DBARRIER_OUTPUT=0

# Debugging
#PARDEBUG = -g3

# Do parallel thetas
#PARTHETAS = -DDO_PARALLELTHETAS=0

# Turn off loadbalancing
# not needed now because set to 0 in commondefs.h
#DO_LOAD = -DDO_LOADBALANCE=0

# Suppress loadbalance output
#LOAD_OUT = -DLOADBALANCE_OUTPUT=0

# Compiler warnings
#WARNFLAG = -Wall


############## parallel CFLAGS ########################

# Required parallel flags
PARFLAGS = \
	$(BAR_OUT) \
	$(PARDEBUG) \
	$(DO_LOAD) \
	$(DOPAR) \
	$(EXTRAFLAGS) \
	$(LOAD_OUT) \
	$(OPTFLAG) \
	$(NOCKPT) \
	$(PARINCLPATH) \
	$(PARLIB) \
	$(PARTHETAS) \
	$(PROTO) \
	$(SYSDEP) \
	$(WARNFLAG)

# parallel ILINK CFLAGS
ILCFLAGS = \
	$(ILPARMEM) \
	$(ILPRECOMP) \
	$(PARFLAGS)

# parallel MLINK CFLAGS
MLCFLAGS = \
	$(MLPARMEM) \
	$(MLPRECOMP) \
	$(PARFLAGS)

# parallel LINKMAP CFLAGS
LICFLAGS = \
	$(LIPARMEM) \
	$(LIPRECOMP) \
	$(PARFLAGS)

# Program-specific CFLAGS
PARILCFLAGS = $(ILCFLAGS) -DILINK
PARMLCFLAGS = $(MLCFLAGS) -DMLINK
PARLICFLAGS = $(LICFLAGS) -DLINKMAP


######################################################
################# common object files ################
######################################################


# Object files common to all sequential programs,
# each of which has a "slow" counterpart.
FASTCOMOBJS = \
	automodified.o \
	comlike.o \
	sexmodified.o

# "Slow" versions of FASTCOMOBJS.
SLOWCOMOBJS = \
	slowautomodified.o \
	slowcomlike.o \
	slowsexmodified.o


######################################################
################# ILINK object files #################
######################################################

# Object files common to "fast" and "slow" versions of ILINK.
ILCOMOBJS = \
	iliostuff.o \
	ilcheck.o \
	ilcominputcode.o \
	ilcomrec.o \
	ilgetvect.o \
	ilinputcode.o \
	iloldnuclear.o

# "Fast" versions of ILINK-specific object files that
# have corresponding "slow" versions
FASTILOBJS = \
	ilcommoncode.o \
	ilink.o

# Object files for "fast" ILINK.
ILOBJS = \
	$(ILCOMOBJS) \
	$(FASTCOMOBJS) \
	$(FASTILOBJS)

# Object files for "slow" ILINK
ILSLOWOBJS = \
	$(ILCOMOBJS) \
	$(SLOWCOMOBJS) \
	slowilcommoncode.o \
	slowilink.o

# Object files for parallel ILINK
PARILOBJS = \
	parilcominputcode.o \
	parilcomlike.o \
	parilcommoncode.o \
	parilcompar.o \
	parilcomrec.o \
	parilgetvect.o \
	parilink.o \
	parilinputcode.o \
	pariliostuff.o \
	parilmodified.o \
	pariloldnuclear.o


######################################################
################ LINKMAP object files ################
######################################################

# Object files common to "fast" and "slow" versions of LINKMAP.
LICOMOBJS = \
	liiostuff.o \
	licheck.o \
	licominputcode.o \
	licomrec.o \
	ligetvect.o \
	liinputcode.o \
	lioldnuclear.o

# "Fast" versions of LINKMAP-specific object files that
# have corresponding "slow" versions
FASTLIOBJS = \
	licommoncode.o \
	linkmap.o

# Object files for "fast" LINKMAP.
LIOBJS = \
	$(LICOMOBJS) \
	$(FASTCOMOBJS) \
	$(FASTLIOBJS)

# Object files for "slow" LINKMAP
LISLOWOBJS = \
	$(LICOMOBJS) \
	$(SLOWCOMOBJS) \
	slowlicommoncode.o \
	slowlinkmap.o

# Object files for parallel LINKMAP
PARLIOBJS = \
	parlicheck.o \
	parlicominputcode.o \
	parlicomlike.o \
	parlicommoncode.o \
	parlicomoutput.o \
	parlicompar.o \
	parlicomrec.o \
	parligetvect.o \
	parliinputcode.o \
	parliiostuff.o \
	parlimodified.o \
	parlinkmap.o \
	parlioldnuclear.o \
	parlioutbuff.o \
	parlistrbuff.o


######################################################
################# MLINK object files #################
######################################################

# Object files common to "fast" and "slow" versions of MLINK.
MLCOMOBJS = \
	mliostuff.o \
	mlcheck.o \
	mlcominputcode.o \
	mlcomrec.o \
	mlgetvect.o \
	mlinputcode.o \
	mloldnuclear.o

# "Fast" versions of MLINK-specific object files that
# have corresponding "slow" versions
FASTMLOBJS = \
	mlcommoncode.o \
	mlink.o

# Object files for "fast" MLINK.
MLOBJS = \
	$(MLCOMOBJS) \
	$(FASTCOMOBJS) \
	$(FASTMLOBJS)

# Object files for "slow" MLINK
MLSLOWOBJS = \
	$(MLCOMOBJS) \
	$(SLOWCOMOBJS) \
	slowmlcommoncode.o \
	slowmlink.o

# Object files for parallel MLINK
PARMLOBJS = \
	parmlcheck.o \
	parmlcominputcode.o \
	parmlcomlike.o \
	parmlcommoncode.o \
	parmlcomoutput.o \
	parmlcompar.o \
	parmlcomrec.o \
	parmlgetvect.o \
	parmlink.o \
	parmlinputcode.o \
	parmliostuff.o \
	parmlmodified.o \
	parmloldnuclear.o \
	parmloutbuff.o \
	parmlstrbuff.o


######################################################
############### LODSCORE object files ################
######################################################

# Object files common to "fast" and "slow" versions of LODSCORE.
LOCOMOBJS = \
	loiostuff.o \
	locheck.o \
	locomrec.o \
	logetvect.o \
	loinputcode.o \
	looldnuclear.o

LOOBJS = \
	$(LOCOMOBJS) \
	automodified.o \
	locomlike.o \
	sexmodified.o \
	lodscore.o

LOSLOWOBJS = \
	$(LOCOMOBJS) \
	slowautomodified.o \
	slowlocomlike.o \
	slowsexmodified.o \
	slowlodscore.o


######################################################
#################### targets #########################
######################################################

# Executables
FASTEXE = lodscore ilink linkmap mlink unknown
SLOWEXE = slowlodscore slowilink slowlinkmap slowmlink unknown
UDPEXE = ilink.udp linkmap.udp mlink.udp unknown
ATMEXE = ilink.atm linkmap.atm mlink.atm unknown
P4EXE = ilink.p4 linkmap.p4 mlink.p4 unknown

CLEANEXE = \
	$(BINDIR)/lodscore \
	$(BINDIR)/ilink \
	$(BINDIR)/linkmap \
	$(BINDIR)/mlink \
	$(BINDIR)/slowlodscore \
	$(BINDIR)/slowilink \
	$(BINDIR)/slowlinkmap \
	$(BINDIR)/slowmlink \
	$(BINDIR)/ofm \
	$(BINDIR)/unknown


installfast: $(FASTEXE)

installslow: $(SLOWEXE)

installudp: $(UDPEXE)

installatm: $(ATMEXE)

installp4: $(P4EXE)


################## sequential targets ###############

ilink: $(ILOBJS) 
	$(CC) $(CFLAGS) $(ILOBJS) -o $(BINDIR)/ilink $(LDFLAGS)

slowilink: $(ILSLOWOBJS) 
	$(CC) $(CFLAGS) $(SFLAG) $(ILSLOWOBJS) -o $(BINDIR)/ilink $(LDFLAGS)


lodscore: $(LOOBJS) 
	$(CC) $(CFLAGS) $(LOOBJS) -o $(BINDIR)/lodscore $(LDFLAGS)

slowlodscore: $(LOSLOWOBJS) 
	$(CC) $(CFLAGS) $(SFLAG) $(LOSLOWOBJS) -o $(BINDIR)/lodscore $(LDFLAGS)


linkmap: $(LIOBJS) 
	$(CC) $(CFLAGS) $(LIOBJS) -o $(BINDIR)/linkmap $(LDFLAGS)

slowlinkmap: $(LISLOWOBJS) 
	$(CC) $(CFLAGS) $(SFLAG) $(LISLOWOBJS) -o $(BINDIR)/linkmap $(LDFLAGS)


mlink: $(MLOBJS) 
	$(CC) $(CFLAGS) $(MLOBJS) -o $(BINDIR)/mlink $(LDFLAGS)

slowmlink: $(MLSLOWOBJS) 
	$(CC) $(CFLAGS) $(SFLAG) $(MLSLOWOBJS) -o $(BINDIR)/mlink $(LDFLAGS)


ofm: ofm.c 
	$(CC) $(SMALLCFLAGS) ofm.c -o $(BINDIR)/ofm

loopbrk.o: loopbrk.c unknown.h
	$(CC) $(SMALLCFLAGS) loopbrk.c -c -o loopbrk.o 
unknown.o: unknown.c unknown.h
	$(CC) $(SMALLCFLAGS) unknown.c  -c -o unknown.o  
unknown: unknown.o loopbrk.o
	$(CC) $(SMALLCFLAGS) loopbrk.o unknown.o -o $(BINDIR)/unknown $(LDFLAGS)


################## parallel targets ###############


ilink.atm: $(PARILOBJS) $(TmkLIB)/libTmk_atm.a 
	$(PARCC) $(PARILCFLAGS) -o $(BINDIR)/$(@) $(PARILOBJS) $(ATMLD)

ilink.udp: $(PARILOBJS) $(TmkLIB)/libTmk_udp.a 
	$(PARCC) $(PARILCFLAGS) -o $(BINDIR)/$(@) $(PARILOBJS) $(UDPLD)

ilink.p4: $(PARILOBJS) $(P4_LIBS) 
	$(PARCC) $(PARILCFLAGS) -o $(BINDIR)/$(@) $(PARILOBJS) $(P4LD)
	

linkmap.atm: $(PARLIOBJS) $(TmkLIB)/libTmk_atm.a 
	$(PARCC) $(PARLICFLAGS) -o $(BINDIR)/$(@) $(PARLIOBJS) $(ATMLD)
	
linkmap.udp: $(PARLIOBJS) $(TmkLIB)/libTmk_udp.a 
	$(PARCC) $(PARLICFLAGS) -o $(BINDIR)/$(@) $(PARLIOBJS) $(UDPLD)
	
linkmap.p4: $(PARLIOBJS) $(P4_LIBS) 
	$(PARCC) $(PARLICFLAGS) -o $(BINDIR)/$(@) $(PARLIOBJS) $(P4LD)
	

mlink.atm: $(PARMLOBJS) $(TmkLIB)/libTmk_atm.a 
	$(PARCC) $(PARMLCFLAGS) -o $(BINDIR)/$(@) $(PARMLOBJS) $(ATMLD)

mlink.udp: $(PARMLOBJS) $(TmkLIB)/libTmk_udp.a 
	$(PARCC) $(PARMLCFLAGS) -o $(BINDIR)/$(@) $(PARMLOBJS) $(UDPLD)

mlink.p4: $(PARMLOBJS) $(P4_LIBS) 
	$(PARCC) $(PARMLCFLAGS) -o $(BINDIR)/$(@) $(PARMLOBJS) $(P4LD)
	


######################################################
############## implicit rules for object files #######
######################################################

############## ILINK rules

ilcheck.o:
	$(CC) $(CFLAGS) -DILINK -c check.c -o ilcheck.o

ilcominputcode.o:
	$(CC) $(CFLAGS) -DILINK -c cominputcode.c -o ilcominputcode.o

ilcommoncode.o:
	$(CC) $(CFLAGS) -DILINK -c commoncode.c -o ilcommoncode.o

ilcomrec.o:
	$(CC) $(CFLAGS) -DILINK -c comrec.c -o ilcomrec.o

ilgetvect.o:
	$(CC) $(CFLAGS) -DILINK -c ilgetvect.c -o ilgetvect.o

ilink.o:
	$(CC) $(CFLAGS) -DILINK -c ilink.c -o ilink.o

ilinputcode.o:
	$(CC) $(CFLAGS) -DILINK -c ilinputcode.c -o ilinputcode.o

iliostuff.o:
	$(CC) $(CFLAGS) -DILINK -c iostuff.c -o iliostuff.o

iloldnuclear.o:
	$(CC) $(CFLAGS) -DILINK -c iloldnuclear.c -o iloldnuclear.o

slowilcommoncode.o:
	$(CC) $(CFLAGS) $(SFLAG) -DILINK -c commoncode.c -o slowilcommoncode.o

slowilink.o:
	$(CC) $(CFLAGS) $(SFLAG) -c ilink.c -o slowilink.o

parilcominputcode.o:
	$(PARCC) $(PARILCFLAGS) -c cominputcode.c -o parilcominputcode.o

parilcomlike.o:
	$(PARCC) $(PARILCFLAGS) -c comlike.c -o parilcomlike.o

parilcommoncode.o:
	$(PARCC) $(PARILCFLAGS) -c commoncode.c -o parilcommoncode.o

parilcompar.o:
	$(PARCC) $(PARILCFLAGS) -c compar.c -o parilcompar.o

parilcomrec.o:
	$(PARCC) $(PARILCFLAGS) -c comrec.c -o parilcomrec.o

parilgetvect.o:
	$(PARCC) $(PARILCFLAGS) -c ilgetvect.c -o parilgetvect.o

parilink.o:
	$(PARCC) $(PARILCFLAGS) -c ilink.c -o parilink.o

parilinputcode.o:
	$(PARCC) $(PARILCFLAGS) -c ilinputcode.c -o parilinputcode.o

pariliostuff.o:
	$(PARCC) $(PARILCFLAGS) -c iostuff.c -o pariliostuff.o

parilmodified.o:
	$(PARCC) $(PARILCFLAGS) -c parmodified.c -o parilmodified.o

pariloldnuclear.o:
	$(PARCC) $(PARILCFLAGS) -c iloldnuclear.c -o pariloldnuclear.o



############## LINKMAP rules


licheck.o:
	$(CC) $(CFLAGS) -DLINKMAP -c check.c -o licheck.o

licominputcode.o:
	$(CC) $(CFLAGS) -DLINKMAP -c cominputcode.c -o licominputcode.o

licommoncode.o:
	$(CC) $(CFLAGS) -DLINKMAP -c commoncode.c -o licommoncode.o

licomrec.o:
	$(CC) $(CFLAGS) -DLINKMAP -c comrec.c -o licomrec.o

ligetvect.o:
	$(CC) $(CFLAGS) -DLINKMAP -c ligetvect.c -o ligetvect.o

liinputcode.o:
	$(CC) $(CFLAGS) -DLINKMAP -c liinputcode.c -o liinputcode.o

liiostuff.o:
	$(CC) $(CFLAGS) -DLINKMAP -c iostuff.c -o liiostuff.o

linkmap.o:
	$(CC) $(CFLAGS) -DLINKMAP -c linkmap.c -o linkmap.o

lioldnuclear.o:
	$(CC) $(CFLAGS) -DLINKMAP -c lioldnuclear.c -o lioldnuclear.o

slowlicommoncode.o:
	$(CC) $(CFLAGS) $(SFLAG) -DLINKMAP -c commoncode.c -o slowlicommoncode.o

slowlinkmap.o:
	$(CC) $(CFLAGS) $(SFLAG) -c linkmap.c -o slowlinkmap.o

parlicheck.o:
	$(PARCC) $(PARLICFLAGS) -c check.c -o parlicheck.o

parlicominputcode.o:
	$(PARCC) $(PARLICFLAGS) -c cominputcode.c -o parlicominputcode.o

parlicomlike.o:
	$(PARCC) $(PARLICFLAGS) -c comlike.c -o parlicomlike.o

parlicommoncode.o:
	$(PARCC) $(PARLICFLAGS) -c commoncode.c -o parlicommoncode.o

parlicomoutput.o:
	$(PARCC) $(PARLICFLAGS) -c comoutput.c -o parlicomoutput.o

parlicompar.o:
	$(PARCC) $(PARLICFLAGS) -c compar.c -o parlicompar.o

parlicomrec.o:
	$(PARCC) $(PARLICFLAGS) -c comrec.c -o parlicomrec.o

parligetvect.o:
	$(PARCC) $(PARLICFLAGS) -c ligetvect.c -o parligetvect.o

parliinputcode.o:
	$(PARCC) $(PARLICFLAGS) -c liinputcode.c -o parliinputcode.o

parliiostuff.o:
	$(PARCC) $(PARLICFLAGS) -c iostuff.c -o parliiostuff.o

parlimodified.o:
	$(PARCC) $(PARLICFLAGS) -c parmodified.c -o parlimodified.o

parlinkmap.o:
	$(PARCC) $(PARLICFLAGS) -c linkmap.c -o parlinkmap.o

parlioldnuclear.o:
	$(PARCC) $(PARLICFLAGS) -c lioldnuclear.c -o parlioldnuclear.o

parlioutbuff.o:
	$(PARCC) $(PARLICFLAGS) -c outbuff.c -o parlioutbuff.o

parlistrbuff.o:
	$(PARCC) $(PARLICFLAGS) -c strbuff.c -o parlistrbuff.o


############## MLINK rules


mlcheck.o:
	$(CC) $(CFLAGS) -DMLINK -c check.c -o mlcheck.o

mlcominputcode.o:
	$(CC) $(CFLAGS) -DMLINK -c cominputcode.c -o mlcominputcode.o

mlcommoncode.o:
	$(CC) $(CFLAGS) -DMLINK -c commoncode.c -o mlcommoncode.o

mlcomrec.o:
	$(CC) $(CFLAGS) -DMLINK -c comrec.c -o mlcomrec.o

mlgetvect.o:
	$(CC) $(CFLAGS) -DMLINK -c mlgetvect.c -o mlgetvect.o

mlink.o:
	$(CC) $(CFLAGS) -DMLINK -c mlink.c -o mlink.o

mlinputcode.o:
	$(CC) $(CFLAGS) -DMLINK -c mlinputcode.c -o mlinputcode.o

mliostuff.o:
	$(CC) $(CFLAGS) -DMLINK -c iostuff.c -o mliostuff.o

mloldnuclear.o:
	$(CC) $(CFLAGS) -DMLINK -c mloldnuclear.c -o mloldnuclear.o

slowmlcommoncode.o:
	$(CC) $(CFLAGS) $(SFLAG) -DMLINK -c commoncode.c -o slowmlcommoncode.o

slowmlink.o:
	$(CC) $(CFLAGS) $(SFLAG) -c mlink.c -o slowmlink.o

parmlcheck.o:
	$(PARCC) $(PARMLCFLAGS) -c check.c -o parmlcheck.o

parmlcominputcode.o:
	$(PARCC) $(PARMLCFLAGS) -c cominputcode.c -o parmlcominputcode.o

parmlcomlike.o:
	$(PARCC) $(PARMLCFLAGS) -c comlike.c -o parmlcomlike.o

parmlcommoncode.o:
	$(PARCC) $(PARMLCFLAGS) -c commoncode.c -o parmlcommoncode.o

parmlcomoutput.o:
	$(PARCC) $(PARMLCFLAGS) -c comoutput.c -o parmlcomoutput.o

parmlcompar.o:
	$(PARCC) $(PARMLCFLAGS) -c compar.c -o parmlcompar.o

parmlcomrec.o:
	$(PARCC) $(PARMLCFLAGS) -c comrec.c -o parmlcomrec.o

parmlgetvect.o:
	$(PARCC) $(PARMLCFLAGS) -c mlgetvect.c -o parmlgetvect.o

parmlink.o:
	$(PARCC) $(PARMLCFLAGS) -c mlink.c -o parmlink.o

parmlinputcode.o:
	$(PARCC) $(PARMLCFLAGS) -c mlinputcode.c -o parmlinputcode.o

parmliostuff.o:
	$(PARCC) $(PARMLCFLAGS) -c iostuff.c -o parmliostuff.o

parmlmodified.o:
	$(PARCC) $(PARMLCFLAGS) -c parmodified.c -o parmlmodified.o

parmloldnuclear.o:
	$(PARCC) $(PARMLCFLAGS) -c mloldnuclear.c -o parmloldnuclear.o

parmloutbuff.o:
	$(PARCC) $(PARMLCFLAGS) -c outbuff.c -o parmloutbuff.o

parmlstrbuff.o:
	$(PARCC) $(PARMLCFLAGS) -c strbuff.c -o parmlstrbuff.o


############# LODSCORE rules


locheck.o:
	$(CC) $(CFLAGS) -DLODSCORE -c check.c -o locheck.o

locomlike.o:
	$(CC) $(CFLAGS) -DLODSCORE -c comlike.c -o locomlike.o

locomrec.o:
	$(CC) $(CFLAGS) -DLODSCORE -c comrec.c -o locomrec.o

lodscore.o:
	$(CC) $(CFLAGS) -DLODSCORE -c lodscore.c -o lodscore.o

logetvect.o:
	$(CC) $(CFLAGS) -DLODSCORE -c logetvect.c -o logetvect.o

loinputcode.o:
	$(CC) $(CFLAGS) -DLODSCORE -c loinputcode.c -o loinputcode.o

loiostuff.o:
	$(CC) $(CFLAGS) -DLODSCORE -c iostuff.c -o loiostuff.o

looldnuclear.o:
	$(CC) $(CFLAGS) -DLODSCORE -c looldnuclear.c -o looldnuclear.o


slowlocomlike.o:
	$(CC) $(CFLAGS) $(SFLAG) -DLODSCORE -c comlike.c -o slowlocomlike.o

slowlodscore.o:
	$(CC) $(CFLAGS) $(SFLAG) -DLODSCORE -c lodscore.c -o slowlodscore.o



############# generic sequential code rules


automodified.o:
	$(CC) $(CFLAGS) -c automodified.c -o automodified.o

comlike.o:
	$(CC) $(CFLAGS) -c comlike.c -o comlike.o

sexmodified.o:
	$(CC) $(CFLAGS) -c sexmodified.c -o sexmodified.o

slowautomodified.o:
	$(CC) $(CFLAGS) $(SFLAG) -c slowautomodified.c -o slowautomodified.o

slowcomlike.o:
	$(CC) $(CFLAGS) $(SFLAG) -c comlike.c -o slowcomlike.o

slowsexmodified.o:
	$(CC) $(CFLAGS) $(SFLAG) -c slowsexmodified.c -o slowsexmodified.o


######################################################
############# cleanup rules ##########################
######################################################


clean:
		$(RM) mlink *.o

parclean:
	$(RM) \
	$(PARILOBJS) $(PARMLOBJS) $(PARLIOBJS) \
	$(BINDIR)/*.atm $(BINDIR)/*.udp $(BINDIR)/*.p4

allclean: clean parclean
		$(RM) $(CLEANEXE)


######################################################
############## dependencies ##########################
######################################################

automodified.o :  automodified.c commondefs.h moddefs.h 
comlike.o :  comlike.c commondefs.h moddefs.h 
ilcheck.o :  check.c commondefs.h checkpointdefs.h gemdefs.h ildefs.h 
ilcominputcode.o :  cominputcode.c commondefs.h checkpointdefs.h gemdefs.h 
ilcommoncode.o :  commoncode.c commondefs.h gemdefs.h moddefs.h checkpointdefs.h \
  ildefs.h 
ilcomrec.o :  comrec.c commondefs.h checkpointdefs.h 
ilgetvect.o :  ilgetvect.c commondefs.h comgetvect.c 
ilink.o :  ilink.c commondefs.h checkpointdefs.h gemdefs.h ildefs.h moddefs.h 
ilinputcode.o :  ilinputcode.c commondefs.h gemdefs.h ildefs.h 
iliostuff.o :  iostuff.c commondefs.h 
iloldnuclear.o :  iloldnuclear.c commondefs.h ildefs.h comnuclear.c oldsegup.c 
licheck.o :  check.c commondefs.h checkpointdefs.h gemdefs.h ildefs.h 
licominputcode.o :  cominputcode.c commondefs.h checkpointdefs.h gemdefs.h 
licommoncode.o :  commoncode.c commondefs.h gemdefs.h moddefs.h checkpointdefs.h \
  lidefs.h 
licomrec.o :  comrec.c commondefs.h checkpointdefs.h 
ligetvect.o :  ligetvect.c commondefs.h comgetvect.c 
liinputcode.o :  liinputcode.c commondefs.h gemdefs.h lidefs.h 
liiostuff.o :  iostuff.c commondefs.h 
linkmap.o :  linkmap.c commondefs.h gemdefs.h lidefs.h checkpointdefs.h moddefs.h 
lioldnuclear.o :  lioldnuclear.c commondefs.h lidefs.h comnuclear.c oldsegup.c 
locheck.o :  check.c commondefs.h checkpointdefs.h gemdefs.h ildefs.h 
locomlike.o :  comlike.c commondefs.h moddefs.h 
locomrec.o :  comrec.c commondefs.h checkpointdefs.h lodefs.h 
lodscore.o :  lodscore.c commondefs.h checkpointdefs.h gemdefs.h lodefs.h moddefs.h 
logetvect.o :  logetvect.c commondefs.h lodefs.h comgetvect.c 
loinputcode.o :  loinputcode.c commondefs.h checkpointdefs.h gemdefs.h lodefs.h 
loiostuff.o :  iostuff.c commondefs.h lodefs.h 
looldnuclear.o :  looldnuclear.c commondefs.h lodefs.h comnuclear.c oldsegup.c 
mlcheck.o :  check.c commondefs.h checkpointdefs.h gemdefs.h ildefs.h 
mlcominputcode.o :  cominputcode.c commondefs.h checkpointdefs.h gemdefs.h 
mlcommoncode.o :  commoncode.c commondefs.h gemdefs.h moddefs.h checkpointdefs.h \
  mldefs.h 
mlcomrec.o :  comrec.c commondefs.h checkpointdefs.h 
mlgetvect.o :  mlgetvect.c commondefs.h comgetvect.c 
mlink.o :  mlink.c commondefs.h gemdefs.h mldefs.h checkpointdefs.h moddefs.h 
mlinputcode.o :  mlinputcode.c commondefs.h gemdefs.h mldefs.h 
mliostuff.o :  iostuff.c commondefs.h 
mloldnuclear.o :  mloldnuclear.c commondefs.h mldefs.h comnuclear.c oldsegup.c 
parilcominputcode.o :  cominputcode.c commondefs.h  strbuff.h gemdefs.h 
parilcomlike.o :  comlike.c commondefs.h  strbuff.h moddefs.h 
parilcommoncode.o :  commoncode.c commondefs.h  strbuff.h gemdefs.h moddefs.h \
  ildefs.h compar.h 
parilcompar.o :  compar.c commondefs.h  strbuff.h gemdefs.h moddefs.h compar.h 
parilcomrec.o :  comrec.c commondefs.h  strbuff.h 
parilgetvect.o :  ilgetvect.c commondefs.h  strbuff.h comgetvect.c 
parilink.o :  ilink.c commondefs.h  strbuff.h gemdefs.h ildefs.h moddefs.h \
  compar.h 
parilinputcode.o :  ilinputcode.c commondefs.h  strbuff.h gemdefs.h ildefs.h 
pariliostuff.o :  iostuff.c commondefs.h  strbuff.h 
parilmodified.o :  parmodified.c commondefs.h  strbuff.h moddefs.h 
pariloldnuclear.o :  iloldnuclear.c commondefs.h  strbuff.h ildefs.h comnuclear.c \
  oldsegup.c 
parlicheck.o :  check.c commondefs.h  strbuff.h checkpointdefs.h gemdefs.h ildefs.h 
parlicominputcode.o :  cominputcode.c commondefs.h  strbuff.h gemdefs.h 
parlicomlike.o :  comlike.c commondefs.h  strbuff.h moddefs.h 
parlicommoncode.o :  commoncode.c commondefs.h  strbuff.h gemdefs.h moddefs.h \
  lidefs.h compar.h 
parlicomoutput.o :  comoutput.c commondefs.h  strbuff.h gemdefs.h moddefs.h lidefs.h 
parlicompar.o :  compar.c commondefs.h  strbuff.h gemdefs.h moddefs.h lidefs.h \
  compar.h 
parlicomrec.o :  comrec.c commondefs.h  strbuff.h 
parligetvect.o :  ligetvect.c commondefs.h  strbuff.h comgetvect.c 
parliinputcode.o :  liinputcode.c commondefs.h  strbuff.h gemdefs.h lidefs.h 
parliiostuff.o :  iostuff.c commondefs.h  strbuff.h 
parlimodified.o :  parmodified.c commondefs.h  strbuff.h moddefs.h 
parlinkmap.o :  linkmap.c commondefs.h  strbuff.h gemdefs.h lidefs.h moddefs.h \
  compar.h 
parlioldnuclear.o :  lioldnuclear.c commondefs.h  strbuff.h lidefs.h comnuclear.c \
  oldsegup.c 
parlioutbuff.o :  outbuff.c commondefs.h  strbuff.h 
parlistrbuff.o :  strbuff.c strbuff.h commondefs.h  
parmlcheck.o :  check.c commondefs.h  strbuff.h checkpointdefs.h gemdefs.h ildefs.h 
parmlcominputcode.o :  cominputcode.c commondefs.h  strbuff.h gemdefs.h 
parmlcomlike.o :  comlike.c commondefs.h  strbuff.h moddefs.h 
parmlcommoncode.o :  commoncode.c commondefs.h  strbuff.h gemdefs.h moddefs.h \
  mldefs.h compar.h 
parmlcomoutput.o :  comoutput.c commondefs.h  strbuff.h gemdefs.h moddefs.h mldefs.h 
parmlcompar.o :  compar.c commondefs.h  strbuff.h gemdefs.h moddefs.h mldefs.h \
  compar.h 
parmlcomrec.o :  comrec.c commondefs.h  strbuff.h 
parmlgetvect.o :  mlgetvect.c commondefs.h  strbuff.h comgetvect.c 
parmlink.o :  mlink.c commondefs.h  strbuff.h gemdefs.h mldefs.h moddefs.h \
  compar.h 
parmlinputcode.o :  mlinputcode.c commondefs.h  strbuff.h gemdefs.h mldefs.h 
parmliostuff.o :  iostuff.c commondefs.h  strbuff.h 
parmlmodified.o :  parmodified.c commondefs.h  strbuff.h moddefs.h 
parmloldnuclear.o :  mloldnuclear.c commondefs.h  strbuff.h mldefs.h comnuclear.c \
  oldsegup.c 
parmloutbuff.o :  outbuff.c commondefs.h  strbuff.h 
parmlstrbuff.o :  strbuff.c strbuff.h commondefs.h  
sexmodified.o :  sexmodified.c commondefs.h moddefs.h 
slowautomodified.o :  slowautomodified.c commondefs.h 
slowcomlike.o :  comlike.c commondefs.h 
slowilcommoncode.o :  commoncode.c commondefs.h gemdefs.h checkpointdefs.h ildefs.h 
slowilink.o :  ilink.c commondefs.h checkpointdefs.h gemdefs.h ildefs.h 
slowlicommoncode.o :  commoncode.c commondefs.h gemdefs.h checkpointdefs.h lidefs.h 
slowlinkmap.o :  linkmap.c commondefs.h gemdefs.h lidefs.h checkpointdefs.h 
slowlocomlike.o :  comlike.c commondefs.h 
slowlodscore.o :  lodscore.c commondefs.h checkpointdefs.h gemdefs.h lodefs.h 
slowmlcommoncode.o :  commoncode.c commondefs.h gemdefs.h checkpointdefs.h mldefs.h 
slowmlink.o :  mlink.c commondefs.h gemdefs.h mldefs.h checkpointdefs.h 
slowsexmodified.o :  slowsexmodified.c commondefs.h 


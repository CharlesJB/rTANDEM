#include <Rcpp.h> // rTANDEM
#include "dataLoader.h" // rTANDEM

// rTANDEM : Declaration initialy in tandem.cpp where move in tandem.h
#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
#include <algorithm>

#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"


/*
 * windows.h and the definition of ProcessThread are necessary for multithreading
 * in the Windows 32 environment. UNIX platforms that use POSIX threads use an
 * alternate version of this file.
 */
#ifdef MSVC
	#include "windows.h"
	DWORD WINAPI ProcessThread(LPVOID);
	DWORD WINAPI RefineThread(LPVOID);
#else
	#include <pthread.h>
	void* ProcessThread(void *pParam);
	void* RefineThread(void *pParam);
#endif

bool lessThanSpec(const mspectrum &_l,const mspectrum &_r);

bool lessThanSpec(const mspectrum &_l,const mspectrum &_r)
{
	return _l.m_dMH < _r.m_dMH;
}
/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
RcppExport SEXP tandem(SEXP param, SEXP peptide, SEXP saps, SEXP mods, SEXP spectrum); // rTANDEM

/*
 Copyright (C) 2003-2004 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

// File version: 2005-01-23

#include "stdafx.h"
#include "msequence.h"
#include "mspectrum.h"
#include "msequtilities.h"
#include "mscore_tandem.h"
// Factory instance, registers itself with the mscoremanager.
static mscorefactory_tandem factory;
	
mscorefactory_tandem::mscorefactory_tandem()
{
	mscoremanager::register_factory("tandem", this);
}

mplugin* mscorefactory_tandem::create_plugin()
{
	return new mscore_tandem();
}

mscore_tandem::mscore_tandem(void)
{
	m_fScale = 4.0;
	m_pFactorial = new double[64];
	double dFac = 1.0;
	m_pFactorial[0] = 1.0;
	long a = 1;
	while(a < 64)	{
		dFac *= (double)a;
		m_pFactorial[a] = dFac;
		a++;
	}
	m_uiSimd = 0;
	m_pfLogs = new float[101];
	float fV = 0.01;
	a = 1;
	m_pfLogs[0] = 0;
	while(a < 101)	{
		m_pfLogs[a] = log(fV);
		fV += 0.01;
		a++;
	}
	m_fLog2 = (float)log(2.0);
	m_fLog10 = (float)(1.0/log(10.0));

#ifdef MSVC
	m_uiSimd = check_simd();
	if(m_uiSimd == 1)	{
		cout << "|";
		cout.flush();
	}
#endif
}

mscore_tandem::~mscore_tandem(void)
{
	if(m_pplType)	{
		size_t a = 0;
		while(a < m_vmiType.size())	{
			if(m_pplType[a])	{
				delete m_pplType[a];
			}
			a++;
		}
		delete m_pplType;
	}
	if(m_pfLogs != NULL)	{
		delete m_pfLogs;
	}
	if(m_pFactorial != NULL)
		delete m_pFactorial;
}

/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_tandem::add_mi(mspectrum &_s)
{
	if (!mscore::add_mi(_s))
		return false;

	if (m_vmiType.size() == 0)	{
		m_vmiType.reserve((long)m_vSpec.size());
		m_pplType = new unsigned long *[m_vSpec.size()+1];
		size_t a = 0;
		while(a < m_vSpec.size()+1)	{
			m_pplType[a] = NULL;
			a++;
		}
	}
/*
 * use blur to improve accuracy at the edges of the fragment ion m/z error range
 */
	blur(_s.m_vMI);

	return true;
}

/*
 * blur takes an SMap object and adds values to it adjacent to the 
 * values in the existing map. this process makes the conversion from
 * floating point m/z values to integer SMap index values more accurate.
 * the half-width of the distribution around each initial value is set
 * by the m_fWidth value. For example, 
 *
 * if there is an intensity value I at index value N, 
 * then with a width of 2, the following new values are created,
 * I(N-2),I(N-1),I(N+1),I(N+2)
 *
 * the value for the width has an impact on the performance of the
 * protein modeling session: the wider the width, the more memory that is required to hold
 * the SMaps for the spectra. as that memory expands, the ability of
 * the various processor caches to hold the SMaps will change. keeping
 * the SMaps in the fastest cache speeds up the calculation considerably.
 */
__inline__ bool mscore_tandem::blur(vector<mi> &_s)
{
	vmiType vType;
	MIType uType;
	uType.m_fI = 0.0;
	uType.m_lM = 0;
	long lValue = 0;
	long a = 0;
	long w = (long)(0.1+m_fWidth);
/*
 * fConvert and fFactor are only used if the m_lErrorType uses ppm for fragment ion mass errors
 * if ppm is used, the width at m/z = 200.0 is taken as the base width for blurring & the
 * width is scaled by the measured m/z. 
 * NOTE: the m_fErr value used in the ppm case is: (the error in ppm) x 200
 */
	float fConvert = m_fErr/m_fWidth;
	const float fFactor = (float)200.0/fConvert;
	const size_t tSize = _s.size();
	size_t tCount = 0;
	vType.reserve(tSize*3);
/*
 * add additional values based on the m_fWidth setting
 */
	while(tCount < tSize)	{
		if(_s[tCount].m_fI > 0.5)	{
			lValue = (long)(_s[tCount].m_fM/fConvert);
			a = -1*w;
			if(m_lErrorType & T_FRAGMENT_PPM)	{
				a = (long)((float)a * (float)lValue/fFactor - 0.5);
				if(a > -1*w)
					a = -1*w;
			}
			while(a <= w)	{
				if(uType.m_lM == lValue + a)	{
					if(uType.m_fI < _s[tCount].m_fI)	{
						vType.back().m_fI = _s[tCount].m_fI;
					}
				}
				else	{
						uType.m_lM = lValue + a;
						uType.m_fI = _s[tCount].m_fI;
						vType.push_back(uType);
				}
				a++;
			}
		}
		tCount++;
	}
	m_vmiType.push_back(vType);
	return true;
}

/*
 * hfactor returns a factor applied to the score to produce the
 * hyper score, given a number of ions matched.
 */
double mscore_tandem::hfactor(long _l) {
	if(_l > 63)	{
		return 	m_pFactorial[63];
	}
	return m_pFactorial[_l];
}
/*
 * hconvert is use to convert a hyper score to the value used in a score
 * histogram, since the actual scoring distribution may not lend itself
 * to statistical analysis through use of a histogram.
 */
float mscore_tandem::hconvert(float _f) {
	if(_f <= 0.0)
		return 0.0;
	return m_fScale*(float)log10(_f);
}

__inline__ float mscore_tandem::log_10(float _f)	{
	int iV;
	double man = frexp(_f,&iV);
	int iM = (int)(man*100.0 + 0.5);
	return (iV*m_fLog2 + m_pfLogs[iM])*m_fLog10;
}

bool mscore_tandem::clear()
{
	m_vmiType.clear();
	return true;
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses) and m_pfSeq (scoring
 * weight factors). 
 * convolution scores are the sum of the products of the spectrum intensities
 * and the scoring weight factors for a spectrum.
 * hyper scores are the product of the products of the spectrum intensities
 * and the scoring weight factors for a spectrum.
 * the number of predicted ions that correspond to non-zero intensities in
 * the spectrum are tallied and returned.
 */

double mscore_tandem::dot(unsigned long *_v)
{
	float fScore = 0.0;
	float fValue0 = 0.0;
	unsigned long a = 0;
	unsigned long lCount = 0;
	long lType = 0;
	size_t b = 0;
	vector<MIType>::iterator itType = m_vmiType[m_lId].begin();
	// tType and tTypeSize were added in 2006.09.01 to correct a problem
	// created by VC++ 2005. This new version uses a strict bounds checking
	// style for STL iterators, that cause a run time error if an iterator
	// longer than .end() is produced by incrementing the iterator.
	size_t tType = 0;
	const size_t tTypeSize = m_vmiType[m_lId].size();
	unsigned long *pType = NULL;
	if(!m_pplType[m_lId])	{
		m_pplType[m_lId] = new unsigned long[tTypeSize];
		pType = m_pplType[m_lId];
		while(tType < tTypeSize)	{
			pType[tType] = itType->m_lM;
			++itType;
			tType++;
		}
	}
	else	{
		pType = m_pplType[m_lId];
	}
	itType = m_vmiType[m_lId].begin();
	tType = 0;
	size_t tStep = (size_t)(0.5 + (double)tTypeSize/(double)m_lCount);
	if(tStep < 1)	{
		tStep = 1;
	}
	unsigned long lSeq = m_plSeq[a];
	while(lSeq != 0 && tType != tTypeSize)	{
		lType = 0;
		if(*pType < lSeq)	{
			lType = 1;
			// Generally lots more spectrum peaks than sequence peaks.  Trying
			// large steps first helps reduce performance hit for this operation.
			// This is were the iterator bounds checking failed in VC++ 2005.
			// By checking the size first, the iterator is not evaluated and
			// does not produce the failure.
			while (tType+tStep < tTypeSize && *(pType+tStep) < lSeq) {
				tType += tStep;
				pType += tStep;
			}
			do {
				tType++;
				pType++;
			} while(tType < tTypeSize && *(pType) < lSeq);
		}
		else if(*pType > lSeq)	{
			do {
				a++;
				lSeq = m_plSeq[a];
			} while(*pType > lSeq && lSeq != 0);
		}
		if(lSeq == 0 || tType == tTypeSize)	{
			break;
		}
		if(*pType == lSeq)	{
//			if((itType+tType)->m_fI > 0.0 && m_pfSeq[a] > 0.0)	{
				m_pafI[lCount] = (itType+tType)->m_fI;
				m_pafSeq[lCount] = m_pfSeq[a];
				lCount++;
//			}
		}
		if(lType)	{
			a++;
			lSeq = m_plSeq[a];
		}
		else	{
			tType++;
			pType++;
		}
	}
	*_v = lCount;
	if(lCount == 0)	{
		return fScore;
	}
#ifdef MSVC
	if(!m_uiSimd)	{
		for(a = 0; a < lCount; a++)	{
			fScore += m_pafI[a] * m_pafSeq[a];
		}
		return fScore;
	}
	if(lCount < 5)	{
		for(a = 0; a < lCount; a++)	{
			fScore += m_pafI[a] * m_pafSeq[a];
		}
		return fScore;
	}
	// align the arrays to 4 float boundary
	for(a = 0; lCount % 4 != 0; a++)	{
		m_pafI[lCount] = 0.0;
		m_pafSeq[lCount] = 0.0;
		lCount++;
	}
	// create compatible pointers
	__m128* pSum = (__m128*) m_pafSum;
	__m128* pI = (__m128*) m_pafI;
	__m128* pS = (__m128*) m_pafSeq;
	__m128 Sum = _mm_set1_ps(0.0);
	__m128 S1 = _mm_set1_ps(0.0);
	unsigned int n = lCount/4;
	// perform intrinsic calls for SIMD registers

	for(a = 0; a < n; a++)	{
		S1 = _mm_mul_ps(*pI,*pS);
		Sum = _mm_add_ps(Sum,S1);
		pI++;
		pS++;
	}
	m_um128.m = Sum;
	fScore = m_um128.f[0] + m_um128.f[1] + m_um128.f[2] + m_um128.f[3];
#else
	for(a = 0; a < lCount; a++)	{
		fScore += m_pafI[a] * m_pafSeq[a];
	}
#endif
	return (fScore);
}
#ifdef MSVC
unsigned int mscore_tandem::check_simd(void)
{
	unsigned int iSimd = 0;
	int CPUInfo[4] = {-1};
    __cpuid(CPUInfo, 0);
	unsigned int nIds = CPUInfo[0];
    for (unsigned int i=0; i<=nIds; ++i)	{
        __cpuid(CPUInfo, i);
        // Interpret CPU feature information
        if  (i == 1)	{
			if(((CPUInfo[3] >> 26) & 0x01) && sizeof(int *) == 8)	{
				return 1;
			}
        }
    }
	return iSimd;
}
#endif

/*
double mscore_tandem::dot(unsigned long *_v)
{
	float fScore = 0.0;
	float fValue0 = 0.0;
	unsigned long a = 0;
	unsigned long lCount = 0;
	long lType = 0;
	size_t b = 0;
	vector<MIType>::iterator itType = m_vmiType[m_lId].begin();
	// tType and tTypeSize were added in 2006.09.01 to correct a problem
	// created by VC++ 2005. This new version uses a strict bounds checking
	// style for STL iterators, that cause a run time error if an iterator
	// longer than .end() is produced by incrementing the iterator.
	register size_t tType = 0;
	size_t tTypeSize = m_vmiType[m_lId].size();
	register unsigned long *pType = NULL;
	if(!m_pplType[m_lId])	{
		m_pplType[m_lId] = new unsigned long[tTypeSize];
		pType = m_pplType[m_lId];
		while(tType < tTypeSize)	{
			pType[tType] = itType->m_lM;
			++itType;
			tType++;
		}
	}
	else	{
		pType = m_pplType[m_lId];
	}
	tType = 0;
	size_t tStep = (size_t)(0.5+(float)tTypeSize/(float)m_lCount);
	if(tStep == 0)	{
		tStep = 1;
	}
	register unsigned long lSeq = m_plSeq[a];
	while(lSeq != 0 && tType != tTypeSize)	{
		lType = 0;
		if(pType[tType] < lSeq)	{
			lType = 1;
			// Generally lots more spectrum peaks than sequence peaks.  Trying
			// large steps first helps reduce performance hit for this operation.
			// This is were the iterator bounds checking failed in VC++ 2005.
			// By checking the size first, the iterator is not evaluated and
			// does not produce the failure.
			while (tType+tStep < tTypeSize && *(pType+tType+tStep) < lSeq) {
				tType += tStep;
			}
			do {
				tType++;
			} while(tType < tTypeSize && *(pType+tType) < lSeq);
		}
		else if(pType[tType] > lSeq)	{
			do {
				a++;
				lSeq = m_plSeq[a];
			} while(pType[tType] > lSeq && lSeq != 0);
		}
		if(lSeq == 0 || tType == tTypeSize)	{
			break;
		}
		if(pType[tType] == lSeq)	{
			fValue0 = m_vmiType[m_lId][tType].m_fI * m_pfSeq[a];
			if(fValue0 > 0.0)	{
				lCount++;
				fScore += fValue0;
			}
		}
		if(lType)	{
			a++;
			lSeq = m_plSeq[a];
		}
		else	{
			tType++;
		}
	}
	*_v = lCount;	
	return (fScore);
}
*/
float mscore_tandem::ion_check(const unsigned long _v,const size_t _s)
{
	unsigned long a = 0;
	unsigned long lCount = 0;
	long lType = 0;
	size_t b = 0;
	vector<MIType>::iterator itType = m_vmiType[_s].begin();
	vector<MIType>::const_iterator itStart = m_vmiType[_s].begin();
	vector<MIType>::const_iterator itEnd = m_vmiType[_s].end();
	// tType and tTypeSize were added in 2006.09.01 to correct a problem
	// created by VC++ 2005. This new version uses a strict bounds checking
	// style for STL iterators, that cause a run time error if an iterator
	// longer than .end() is produced by incrementing the iterator.
	size_t tTypeSize = m_vmiType[_s].size();
	size_t tType = tTypeSize/2;
	float fV = (float)1.0;
	itType += tType;

	if(itType->m_lM == _v)	{
		fV = itType->m_fI;
		return fV;
	}
	else if(_v > itType->m_lM)	{
		itType++;
		while(itType != itEnd)	{
			if(itType->m_lM == _v)	{
				fV = itType->m_fI;
				return fV;
			}
			if(_v < itType->m_lM)	{
				return fV;
			}
			itType++;
		}
		return fV;
	}
	else	{
		itType--;
		while(itType != itStart)	{
			if(itType->m_lM == _v)	{
				fV = itType->m_fI;
				return fV;
			}
			if(_v > itType->m_lM)	{
				return fV;
			}
			itType--;
		}
		return fV;
	}
	return fV;
}


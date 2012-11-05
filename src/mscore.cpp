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

// File version: 2003-08-01
// File version: 2004-02-01
// File version: 2004-03-01
// File version: 2004-09-01
// File version: 2004-09-30
// File version: 2005-01-01

/*
 * mscore is the class that contains most of the logic for comparing one sequence with
 * many tandem mass spectra. mprocess contains a comprehensive example of how to use an mscore
 * class. mscore has been optimized for speed, without any specific modifications to take
 * advantage of processor architectures.
 */

#include "stdafx.h"
#include <cmath>
#include <float.h>
#include <algorithm>
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"

/*
 * global less than operator for mspectrumdetails classes: to be used in sort operations
 */
bool lessThanDetails(const mspectrumdetails &_l,const mspectrumdetails &_r)
{
	return _l.m_fL < _r.m_fL; 
}
/*
 * global less than operator for mi classes: to be used in sort operations to achieve
 * list ordered from most intense to least intense
 */
bool lessThanMI(const mi &_l,const mi &_r)
{
	return _l.m_fI > _r.m_fI;
}

mscore::mscore(void) :
	m_seqUtil(masscalc::monoisotopic),
	m_seqUtilAvg(masscalc::average)
{
	m_pSeqUtilFrag = &m_seqUtil;	// default to monoisotopic masses for fragment ions

	m_bIsC = false;
	m_bIsN = false;
	m_pSeq = NULL;

	m_lType = T_Y|T_B;
	m_lErrorType = T_PARENT_DALTONS|T_FRAGMENT_DALTONS;
	m_fParentErrPlus = 2.0;
	m_fParentErrMinus = 2.0;
	m_fErr = (float)0.45;
	
	m_fWidth = 1.0;
	m_lMaxCharge = 100;
	m_dSeqMH = -1.0;
	m_lSize = 256;
	m_pfSeq = new float[m_lSize];
	m_plSeq = new unsigned long[m_lSize];
	m_pSeq = new char[m_lSize];
	m_bIsotopeError = false;
	long a = 0;
	while(a < 20)	{
		m_plCount[a] = 0;
		m_pfScore[a] = 0;
		a++;
	}
	m_fMinMass = 0.0;
	m_fMaxMass = 1.0;
	m_bUsePam = false;
	m_fHomoError = 4.5;
	m_dScale = 1.0;
	m_bMini = false;
	m_iCharge = 1;
	m_bPhosphoBias = true;
	m_plAA = new unsigned long[256];
	memset((void *)m_plAA,0,256*sizeof(unsigned long));
	m_lCount = 0;
}

mscore::~mscore(void)
{
	if(m_pfSeq != NULL)
		delete m_pfSeq;
	if(m_plSeq != NULL)
		delete m_plSeq;
	if(m_pSeq != NULL)
		delete m_pSeq;
	if(m_plAA != NULL)
		delete m_plAA;
}

bool mscore::set_mini(const bool _b)
{
	m_bMini = _b;
	return m_bMini;
}

/*
 * allows score object to issue warnings, or set variable based on xml.
 * default implementation does nothing.
 */
bool mscore::load_param(XmlParameter &_x)
{
	string strKey = "spectrum, fragment mass type";
	string strValue;
	_x.get(strKey,strValue);
	if (strValue == "average") {
		set_fragment_masstype(masscalc::average);
	}

	return true;
}

/*
 * called before spectrum conditioning to allow the score object to
 * modify the spectrum in ways specific to the scoring algorithm.
 * default implementation does nothing.
 */
bool mscore::precondition(mspectrum &_s)
{
	return true;
}

/*
 * called before scoring inside the score() function to allow any
 * necessary resetting of member variables.
 */
__inline__ void mscore::prescore(const size_t _i)
{
	m_lId = _i;
	m_fHyper = (float)0.0;
}

bool mscore::clear()
{
	m_vSpec.clear(); // vector of all spectra being considered
	                        // for all spectra being considered
	m_vDetails.clear(); // vector of mspectrumdetails objects, for looking up parent ion M+H
	                                     // values of mass spectra
	return true;
}
/*
 * create list of non-zero predicted intensity values for a-ions and their
 * integer converted m/z values
 */
bool mscore::add_A(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
/*
 * get the conversion factor between a straight sequence mass and an a-ion
 */
	double dValue = m_pSeqUtilFrag->m_dA;
/*
 * deal with protein N-terminus
 */
	if(m_bIsN)	{
		dValue += m_pSeqUtilFrag->m_fNT;		
	}
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	size_t tC = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfAScore;
	unsigned long lCount = 0;
/*
 * from N- to C-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	m_dWE = m_fWidth/m_fErr;
	const double dZ = (double)_c;
	while(a < m_lSeqLength)	{
		tC = m_pSeq[a];
#ifdef PLUGGABLE_SCORING
		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
		lValue = mconvert(dValue, _c);
#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
		lValue = mconvert(dValue, dZ);
#endif
		m_plSeq[lCount] = lValue;
		m_pfSeq[lCount] = pfScore[tC];
		lCount++;
		a++;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
	m_lCount = lCount;
	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for b-ions and their
 * integer converted m/z values
 */
bool mscore::add_B(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
/*
 * get the conversion factor between a straight sequence mass and a b-ion
 */
	double dValue = m_pSeqUtilFrag->m_dB;
/*
 * deal with protein N-terminus
 */
	if(m_bIsN)	{
		dValue += m_pSeqUtilFrag->m_fNT;		
	}
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	long lCount = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfBScore;
	float *pfScorePlus = m_pSeqUtilFrag->m_pfYScore;
/*
 * from N- to C-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfBScore
 */
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	size_t tC = 0;
	m_dWE = m_fWidth/m_fErr;
	const double dZ = (double)_c;
	while(a < m_lSeqLength-1)	{
		tC = m_pSeq[a];
#ifdef PLUGGABLE_SCORING
		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
		lValue = mconvert(dValue, _c);
#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
		lValue = mconvert(dValue, dZ);
#endif
		m_plSeq[lCount] = lValue;
		m_pfSeq[lCount] = pfScore[tC]*pfScorePlus[m_pSeq[a+1]];
		if(a == 1)	{
			if(m_pSeq[1] == 'P')	{
				m_pfSeq[lCount] *= 10;
			}
			else	{
				m_pfSeq[lCount] *= 3;
			}
		}
		lCount++;
		a++;
	}
	m_lCount = lCount;
	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for c-ions and their
 * integer converted m/z values
 */
bool mscore::add_C(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
/*
 * get the conversion factor between a straight sequence mass and a b-ion
 */
	double dValue = m_pSeqUtilFrag->m_dC;
/*
 * deal with protein N-terminus
 */
	if(m_bIsN)	{
		dValue += m_pSeqUtilFrag->m_fNT;		
	}
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	size_t tC = 0;
	long lCount = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfBScore;
	float *pfScorePlus = m_pSeqUtilFrag->m_pfYScore;
/*
 * from N- to C-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfBScore
 */
	m_dWE = m_fWidth/m_fErr;
	double dZ = (double)_c;
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a < m_lSeqLength-2)	{
		tC = m_pSeq[a];
#ifdef PLUGGABLE_SCORING
		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
		lValue = mconvert(dValue, _c);
#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
		lValue = mconvert(dValue, dZ);
#endif
		m_plSeq[lCount] = lValue;
		m_pfSeq[lCount] = pfScore[tC]*pfScorePlus[m_pSeq[a+1]];
		lCount++;
		a++;
	}
	m_lCount = lCount;
	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for x-ions and their
 * integer converted m/z values
 */
bool mscore::add_X(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;
/*
 * get the conversion factor between a straight sequence mass and an x-ion
 */
	double dValue = m_pSeqUtilFrag->m_dX;
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];
/*
 * deal with protein C-teminus
 */
	if(m_bIsC)	{
		dValue += m_pSeqUtilFrag->m_fCT;		
	}
	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	size_t tC = 0;
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfXScore;
/*
 * from C- to N-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	m_dWE = m_fWidth/m_fErr;
	double dZ = (double)_c;
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		tC = m_pSeq[a];
#ifdef PLUGGABLE_SCORING
		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
		lValue = mconvert(dValue, _c);
#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
		lValue = mconvert(dValue, dZ);
#endif
		m_plSeq[lCount] = lValue;
		m_pfSeq[lCount] = pfScore[tC];
		lCount++;
		a--;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
	m_lCount = lCount;
	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for y-ions and their
 * integer converted m/z values
 */
bool mscore::add_Y(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;
/*
 * get the conversion factor between a straight sequence mass and a y-ion
 */
	double dValue = m_pSeqUtilFrag->m_dY;
	unsigned long lValue = 0;
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];
/*
/*
 * deal with protein C-teminus
 */
	if(m_bIsC)	{
		dValue +=  m_pSeqUtilFrag->m_fCT;		
	}
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfYScore;
	float *pfScoreMinus = m_pSeqUtilFrag->m_pfBScore;
/*
 * from C- to N-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	long tPos = (unsigned long) m_tSeqPos;
	size_t tC = 0;
	m_dWE = m_fWidth/m_fErr;
	double dZ = (double)_c;
	bool bZero = false;
	if(_t == 0)	{
		bZero = true;
	}
	while(a > 0)	{
		tC = m_pSeq[a];
#ifdef PLUGGABLE_SCORING
		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
		lValue = mconvert(dValue, _c);
#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
		lValue = mconvert(dValue, dZ);
#endif
		if(bZero)	{
			if(a < 5)	{
				m_plSeq[lCount] = lValue;
				m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
				lCount++;
			}
		}
		else	{
			m_plSeq[lCount] = lValue;
			m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
			if(a == 2)	{
				if(m_pSeq[1] == 'P')	{
					m_pfSeq[lCount] *= 10;
				}
				else	{
					m_pfSeq[lCount] *= 3;
				}
			}
			lCount++;
		}
		a--;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
	m_lCount = lCount;
	m_plSeq[lCount] = 0;
	return true;
}

/*
 * create list of non-zero predicted intensity values for y-ions and their
 * integer converted m/z values
 */
bool mscore::add_Z(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;
/*
 * get the conversion factor between a straight sequence mass and a y-ion
 */
	double dValue = m_pSeqUtilFrag->m_dZ;
	unsigned long lValue = 0;
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];
/*
/*
 * deal with protein C-teminus
 */
	if(m_bIsC)	{
		dValue +=  m_pSeqUtilFrag->m_fCT;		
	}
	size_t tC = 0;
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfYScore;
	float *pfScoreMinus = m_pSeqUtilFrag->m_pfBScore;
/*
 * from C- to N-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	m_dWE = m_fWidth/m_fErr;
	double dZ = (double)_c;
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		tC = m_pSeq[a];
#ifdef PLUGGABLE_SCORING
		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
		lValue = mconvert(dValue, _c);
#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
		lValue = mconvert(dValue, dZ);
#endif
		m_plSeq[lCount] = lValue;
		m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
		lCount++;
		m_plSeq[lCount] = mconvert(dValue+m_pSeqUtilFrag->m_dHydrogen, dZ);
		m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
		lCount++;
		a--;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
	m_lCount = lCount;
	m_plSeq[lCount] = 0;
	return true;
}

bool mscore::sort_details()
{
/*
 * update the mstate object
 */
	m_State.create_equals((long)m_vSpec.size());
/*
 * sort the m_vDetails vector to improve efficiency at modeling the vector
 */
	sort(m_vDetails.begin(),m_vDetails.end(),lessThanDetails);
/*
 * store the mspectrumdetails object for the mspectrum
 */

	m_lDetails = (long)m_vDetails.size();
	size_t tLimit = m_vDetails.size();
	size_t a = 0;
	m_sIndex.clear();
	mspectrumindex indTemp;
	float fLast = 0.0;
	while(a < tLimit)	{
		indTemp.m_fM = m_vDetails[a].m_fU;
		indTemp.m_tA = (unsigned long)a;
		if(indTemp.m_fM != fLast)	{
			m_sIndex.insert(indTemp);
		}
		fLast = indTemp.m_fM;
		a += 5;
	}
	return true;
}
/*
 * add_mi does the work necessary to set up an mspectrum object for modeling.
 * default implementation simply checks for errors, and sets spectrum count.
 * override in an mscore derived class implementation to do algorithm specific
 * processing.
 */
bool mscore::add_mi(mspectrum &_s)
{
/*
 * the fragment ion error cannot be zero
 */
	if(m_fErr == 0.0)
		return false;

	m_lSpectra = (long)m_vSpec.size();
	return true;
}

/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - a copy of the mspectrum object is added to m_vSpec
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore::add_details(mspectrum &_s)
{
/*
 * if there is a limit on the number of peaks, sort the temporary mspectrum.m_vMI member
 * using lessThanMI (see top of this file). the sort results in m_vMI having the most
 * intense peaks first in the vector. then, simply erase all but the top m_lMaxPeaks values
 * from that vector and continue.
/*
 * the fragment ion error cannot be zero
 */
	if(m_fErr == 0.0)
		return false;
/*
 * create a temporary mspec object
 */
	mspec spCurrent;
	spCurrent = _s;
/*
 * store the mspec object
 */
	m_vSpec.push_back(spCurrent);
	mspectrumdetails detTemp;
	if(m_lErrorType & T_PARENT_PPM)	{
		detTemp.m_fL = (float)(_s.m_dMH - (_s.m_dMH*m_fParentErrPlus/1e6));
		detTemp.m_fU = (float)(_s.m_dMH + (_s.m_dMH*m_fParentErrMinus/1e6));
	}
	else	{
		detTemp.m_fL = (float)(_s.m_dMH - m_fParentErrPlus);
		detTemp.m_fU = (float)(_s.m_dMH + m_fParentErrMinus);
	}
	if(detTemp.m_fU > m_fMaxMass)	{
		m_fMaxMass = detTemp.m_fU;
	}
	detTemp.m_lA = (unsigned long)m_vSpec.size() - 1;
//	if(m_bIsotopeError)	{
//			detTemp.m_fL += (float)(1.00335);
//			detTemp.m_fU += (float)(1.00335);
//			m_vDetails.push_back(detTemp);
//			detTemp.m_fL -= (float)(1.00335);
//			detTemp.m_fU -= (float)(1.00335);
//	}
	m_vDetails.push_back(detTemp);
	if(m_bIsotopeError)	{
		if(spCurrent.m_fMH > 1000.0)	{
			detTemp.m_fL -= (float)(1.008664916);
			detTemp.m_fU -= (float)(1.008664916);
			m_vDetails.push_back(detTemp);
		}
		if(spCurrent.m_fMH > 1500.0)	{
			detTemp.m_fL -= (float)1.008664916;
			detTemp.m_fU -= (float)1.008664916;
			m_vDetails.push_back(detTemp);
		}
	}
	return true;
}
/*
 * add_seq stores a sequence, if the current value of m_pSeq contains the N-terminal portion
 * of the new sequence. this method is part of optimizing the scoring process: all references
 * to it could be replaced by set_seq. set_seq must be called before add_seq is called.
 */
unsigned long mscore::add_seq(const char *_s,const bool _n,const bool _c,const unsigned long _l,const int _f)
{
	m_tSeqPos = _f;
	if(_s == NULL)
		return 0;
	unsigned long lOldLength = m_lSeqLength;
	m_lSeqLength = _l;
/*
 * if the sequence is too long, use set_seq to adjust the arrays and store the sequence
 */
	if(m_lSeqLength >= m_lSize-1)	{
		return set_seq(_s,_n,_c,_l,_f);
	}
/*
 * copy the new part of the sequence into m_pSeq 
 */
	strcpy(m_pSeq + lOldLength,_s + lOldLength);
	unsigned long a = lOldLength;
	m_bIsC = _c;
	m_State.initialize(m_pSeq,m_lSize);
	m_Term.initialize(m_seqUtil.m_pdAaMod['['],m_seqUtil.m_pdAaMod[']']);
/*
 * update the parent ion M+H value
 */
	map<size_t,size_t>::iterator itValue;
	map<size_t,size_t>::iterator itEnd = m_seqUtil.m_mapMotifs.end();
	SMap::iterator itSeq;
	SMap::iterator itSeqEnd = m_seqUtil.m_mapMods.end();
	if(m_seqUtil.m_bPotentialMotif)	{
		m_seqUtil.clear_motifs(false);
	}
	while(a < m_lSeqLength)	{
		m_dSeqMH += m_seqUtil.m_pdAaMass[m_pSeq[a]] + m_seqUtil.m_pdAaMod[m_pSeq[a]] + m_seqUtil.m_pdAaFullMod[m_pSeq[a]];
		if(m_seqUtil.m_bSequenceMods)	{
			itSeq = m_seqUtil.m_mapMods.find((unsigned long)(m_tSeqPos+a));
			if(itSeq != itSeqEnd)	{
				m_dSeqMH += itSeq->second;
			}
		}
		if(m_seqUtil.m_pdAaMod[m_pSeq[a]+32] != 0.0)	{
			m_State.add_mod(m_pSeq+a);
		}
		if(m_seqUtil.m_bPotentialMotif)	{
			itValue = m_seqUtil.m_mapMotifs.find(m_tSeqPos+a);
			if(itValue != itEnd){
				m_State.add_mod(m_pSeq+a);
				m_seqUtil.add_mod(m_pSeq[a],itValue->second);
			}
		}
		a++;
	}
	if(m_seqUtil.m_bPotentialMotif)	{
		m_seqUtil.set_motifs();
	}
/*
 * deal with protein terminii
 */
	if(m_bIsC)
		m_dSeqMH += m_seqUtil.m_fCT;
/*
 * update the mstate object 
 */
	m_State.m_dSeqMHS = m_dSeqMH;
	m_fMinMass = (float)m_dSeqMH;
	if(m_bUsePam)	{
		m_Pam.initialize(m_pSeq,(size_t)m_lSize,(float)m_dSeqMH);
	}
	if(m_bUseSaps)	{
		m_Sap.initialize(m_pSeq,(size_t)m_lSize,(float)m_dSeqMH);
	}
	return m_lSeqLength;
}
/*
 * check_parents is used by the state machine for dealing with potentially
 * modified sequences. it determines how many parent ions are eligible
 * to be generated by the current modified sequence and stores the
 * vector indices of those eligible spectra in the m_State object.
 * this test is done to improve performance. inlining is not necessary,
 * but this method is called very often in a normal protein modeling session, so removing
 * the method calling overhead can improve overall performance
 */
__inline__ bool mscore::check_parents(void)	{
/*
 * this check improves performance because of the way the state machine assigns modifications
 * if the state machine sequence order is modified, this mechanism should be reviewed
 */
	if(m_State.m_dSeqMHFailedS == m_dSeqMH)	{
		m_State.m_lEqualsS = 0;
		return false;
	}
	if(m_dSeqMH < m_vDetails[0].m_fL)	{
		return false;
	}
	if(m_dSeqMH > m_vDetails[m_lDetails-1].m_fU)	{
		return false;
	}
	vector<mspectrumdetails>::iterator itDetails = m_vDetails.begin();
	vector<mspectrumdetails>::iterator itEnd = m_vDetails.end();
	float fSeqMH = (float)m_dSeqMH;
	set<mspectrumindex>::iterator itIndex;
	mspectrumindex indTemp;
	indTemp.m_fM = fSeqMH;
	if(!m_sIndex.empty())	{
		itIndex = m_sIndex.lower_bound(indTemp);
		if(itIndex != m_sIndex.begin())	{
			itIndex--;
		}
		itDetails = itDetails + (*itIndex).m_tA;
	}
	while(itDetails != itEnd)	{
		if(*itDetails == fSeqMH)	{
			m_State.m_lEqualsS = 0;
			m_State.m_plEqualsS[m_State.m_lEqualsS] = itDetails->m_lA;
			m_State.m_lEqualsS++;
			itDetails++;
			while(itDetails != itEnd && *itDetails == fSeqMH)	{
				m_State.m_plEqualsS[m_State.m_lEqualsS] = itDetails->m_lA;
				m_State.m_lEqualsS++;
				itDetails++;
			}
			return true;
		}
		if(fSeqMH < itDetails->m_fL)	{
			break;
		}
		itDetails++;
	}

/*
 * this check improves performance because of the way the state machine assigns modifications
 * if the state machine sequence order is modified, this mechanism should be reviewed
 */
	m_State.m_dSeqMHFailedS = m_dSeqMH;
	m_State.m_lEqualsS = 0;
	return false;
}
/*
 * get_aa is used to check the current sequence and determine how many of 
 * the residues correspond to modification sites. this method is not used
 * by mscore: check mprocess to see an example of how it is used to create
 * a list of modified residues in a domain
 */
bool mscore::get_aa(vector<maa> &_m,const size_t _a,double &_d)
{
	_d = 1000000.0;
	_m.clear();
	size_t a = 0;
	char cRes = '\0';
	char *pValue = NULL;
	maa aaValue;
	double dDelta = 0.0;
	while(a < 128)	{
		if(a == '[' && m_seqUtil.m_pdAaFullMod['['] != 0.0)	{
			cRes = m_pSeq[0];
			aaValue.m_cRes = cRes;
			aaValue.m_dMod = (float)m_seqUtil.m_pdAaFullMod['['];
			aaValue.m_lPos = (unsigned int)_a;
			_m.push_back(aaValue);
		}
		if(a == ']' && m_seqUtil.m_pdAaFullMod[']'] != 0.0)	{
			cRes = m_pSeq[strlen(m_pSeq) - 1];
			aaValue.m_cRes = cRes;
			aaValue.m_dMod = (float)m_seqUtil.m_pdAaFullMod[']'];
			aaValue.m_lPos = (unsigned int)(_a + strlen(m_pSeq) - 1);
			_m.push_back(aaValue);
		}
		if(m_seqUtil.m_pdAaMod[a] != 0.0)	{
			cRes = (char)a;
			pValue = strchr(m_pSeq,cRes);
			while(pValue != NULL)	{
				aaValue.m_cRes = cRes;
				aaValue.m_dMod = (float)m_seqUtil.m_pdAaMod[a];
				aaValue.m_dPrompt = (float)m_seqUtil.m_pdAaPrompt[a];
				aaValue.m_lPos = (unsigned int)(_a + (pValue - m_pSeq));
				dDelta += aaValue.m_dMod;
				_m.push_back(aaValue);
				pValue++;
				pValue = strchr(pValue,cRes);
			}
		}
		if(m_seqUtil.m_pdAaFullMod[a] != 0.0)	{
			cRes = (char)a;
			pValue = strchr(m_pSeq,cRes);
			while(pValue != NULL)	{
				aaValue.m_cRes = cRes;
				aaValue.m_dMod = (float)m_seqUtil.m_pdAaFullMod[a];
				aaValue.m_lPos = (unsigned int)(_a + (pValue - m_pSeq));
				aaValue.m_dPrompt = 0.0;
				_m.push_back(aaValue);
				pValue++;
				pValue = strchr(pValue,cRes);
			}
		}
		a++;
	}		
	if(m_seqUtil.m_bSequenceMods)	{
		SMap::iterator itValue = m_seqUtil.m_mapMods.begin();
		SMap::const_iterator itEnd = m_seqUtil.m_mapMods.end();
		size_t tValue = 0;
		size_t tEnd = m_tSeqPos + m_lSeqLength;
		while(itValue != itEnd)	{
			tValue = itValue->first;
			if(tValue >= m_tSeqPos || tValue < tEnd)	{
				tValue = tValue - m_tSeqPos;
				cRes = m_pSeq[tValue];
				aaValue.m_cRes = cRes;
				aaValue.m_dMod = (float)itValue->second;
				aaValue.m_lPos = (unsigned int)(_a+tValue);
				if(cRes <= 'Z')	{
					cRes += 32;
				}
				dDelta += aaValue.m_dMod;
				aaValue.m_dPrompt = (float)m_seqUtil.m_pdAaPrompt[cRes];
				_m.push_back(aaValue);
			}
			itValue++;
		}
	}
	if(m_Term.m_lN)	{
		aaValue.m_dMod = (float)m_seqUtil.m_pdAaMod['['];
		aaValue.m_dPrompt = 0.0;
		aaValue.m_lPos = (unsigned int)_a;
		aaValue.m_cRes = m_pSeq[0];
		dDelta += aaValue.m_dMod;
		_m.push_back(aaValue);
	}
	if(m_Term.m_lC)	{
		aaValue.m_dMod = (float)m_seqUtil.m_pdAaMod[']'];
		aaValue.m_dPrompt = 0.0;
		aaValue.m_lPos = (unsigned int)(_a+ m_lSeqLength - 1);
		aaValue.m_cRes = m_pSeq[m_lSeqLength - 1];
		dDelta += aaValue.m_dMod;
		_m.push_back(aaValue);
	}
	if(m_Pam.m_tCount > 0)	{
		aaValue.m_dMod = (float)(m_seqUtil.m_pdAaMass[m_pSeq[m_Pam.m_tPos]] - m_seqUtil.m_pdAaMass[m_Pam.m_pSeqTrue[m_Pam.m_tPos]]);
		aaValue.m_dPrompt = (float)(m_seqUtil.m_pdAaPrompt[m_pSeq[m_Pam.m_tPos]] - m_seqUtil.m_pdAaPrompt[m_Pam.m_pSeqTrue[m_Pam.m_tPos]]);
		aaValue.m_lPos = (unsigned int)(_a+ m_Pam.m_tPos);
		aaValue.m_cRes = m_Pam.m_pSeqTrue[m_Pam.m_tPos];
		aaValue.m_cMut = m_pSeq[m_Pam.m_tPos];
		aaValue.m_strId.clear();
		dDelta += aaValue.m_dMod;
		_d = dDelta;
		_m.push_back(aaValue);
	}
	if(m_Sap.m_tCount > 0)	{
		aaValue.m_dMod = (float)(m_seqUtil.m_pdAaMass[m_pSeq[m_Sap.m_tPos]] - m_seqUtil.m_pdAaMass[m_Sap.m_pSeqTrue[m_Sap.m_tPos]]);
		aaValue.m_dPrompt = (float)(m_seqUtil.m_pdAaPrompt[m_pSeq[m_Sap.m_tPos]] - m_seqUtil.m_pdAaPrompt[m_Sap.m_pSeqTrue[m_Sap.m_tPos]]);
		aaValue.m_lPos = (unsigned int)(_a+ m_Sap.m_tPos);
		aaValue.m_cRes = m_Sap.m_pSeqTrue[m_Sap.m_tPos];
		aaValue.m_cMut = m_pSeq[m_Sap.m_tPos];
		aaValue.m_strId = m_Sap.m_strId;
		dDelta += aaValue.m_dMod;
		_d = dDelta;
		_m.push_back(aaValue);
	}
	size_t b = 0;
	while(b < _m.size())	{
		if(_m[b].m_cRes >= 'a' && _m[b].m_cRes <= 'z')	{
			_m[b].m_cRes -= 32;
		}
		if(_m[b].m_cMut >= 'a' && _m[b].m_cMut <= 'z')	{
			_m[b].m_cMut -= 32;
		}
		b++;
	}
	return true;
}
/*
 * mconvert converts from mass and charge to integer ion value
 * for mi vector.
 */
__inline__ unsigned long mscore::mconvert(double _m, const long _c)
{
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	const double dZ = (double)_c;
	return (unsigned long)((m_pSeqUtilFrag->m_dProton + _m/dZ)*m_fWidth/m_fErr);
}
#ifndef PLUGGABLE_SCORING

__inline__ unsigned long mscore::mconvert(double _m, const double _z)
{
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	return (unsigned long)((m_pSeqUtilFrag->m_dProton + _m/_z)*m_dWE);
}
#endif

/*
 * hfactor returns a factor applied to the score to produce the
 * hyper score, given a number of ions matched.
 */
double mscore::hfactor(long _l)
{
	return 1.0;
}
/*
 * sfactor returns a factor applied to the final convolution score.
 */
double mscore::sfactor()
{
	return 1.0;
}
/*
 * hconvert is use to convert a hyper score to the value used in a score
 * histogram, since the actual scoring distribution may not lend itself
 * to statistical analysis through use of a histogram.
 */
float mscore::hconvert(float _f)
{
	if(_f <= 0.0)
		return 0.0;
	return (float)(m_dScale*_f);
}
/*
 * report_score formats a hyper score for output.
 */
void mscore::report_score(char* buffer, float hyper)
{
	sprintf(buffer,"%.1f",hconvert(hyper));
}

/*
 * load_next is used access the potential sequence modification state machine.
 * it runs the state machine until either:
 *   - a modified sequence that has a mass that is within error of 
 *     one of the mspectrum objects has been found, or
 *   - the last of the modified sequences has been returned
 * it returns true if there are more sequences to consider.
 * NOTE: the last sequence generated by the state machine is the one without any
 *       of the potential modifications in place. maintaining this order is important
 *       to the correct functioning of the logic using in mprocess.
 */
bool mscore::load_next(void)
{
	bool bReturn = false;
	if(!(m_bUsePam || m_bUseSaps))	{
		bReturn = load_state();
		if(bReturn)	{
			return bReturn;
		}
		bReturn = load_next_term();
		return bReturn;
	}
	else if(m_bUsePam)	{
// Shift mutation to be checked
		if(m_Pam.m_tCount != 0)	{
			bReturn = load_state();
			if(bReturn)	{
				return bReturn;
			}
			bReturn = load_next_term();
			if(bReturn)	{
				return bReturn;
			}
		}
		bReturn = load_next_pam();
		return bReturn;
	}
	else if(m_bUseSaps)	{
		bReturn = load_state();
		if(bReturn)	{
			return bReturn;
		}
		bReturn = load_next_term();
		if(bReturn)	{
			return bReturn;
		}
		if(m_Sap.m_bOk)	{
			bReturn = load_next_sap();
		}
		return bReturn;
	}
	return bReturn;
}
/*
 * as of version 2004.03.01, a new state machine has been added that allows for testing
 * for N- and C- terminal modifications as partial modifications. Previous versions
 * tested for these modifications as complete modifications only. This state machine
 * relies on the mscoreterm class, m_Term, to maintain information about the current modification
 * state of the N- and C- terminus of the peptide. See the source for that class (in mscorestate.h)
 * for an explanation of how the state is stored and referenced.
 */
bool mscore::load_next_term(void)
{
	if(!m_Term.m_bC && !m_Term.m_bN)	{
		return false;
	}
// reset the state machine, as there are no further states
	if(m_Term.m_lState == 3)	{
		m_Term.initialize(m_seqUtil.m_pdAaMod['['],m_seqUtil.m_pdAaMod[']']);
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMod['['];
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMod[']'];
		m_State.initialize(m_pSeq,m_lSeqLength);
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float)m_dSeqMH;
		m_State.m_lEqualsS = 0;
		check_parents();
		return false;
	}
	bool bReturn = false;
	if(m_Term.m_lState == 0)	{
// deal with modification at N-terminus, if possible
		if(m_Term.m_bN)	{
			m_Term.m_lState = 1;
			m_dSeqMH += m_pSeqUtilFrag->m_pdAaMod['['];
			m_Term.m_lN = 1;
			m_Term.m_lC = 0;
			m_State.initialize(m_pSeq,m_lSeqLength);
			m_State.m_dSeqMHS = m_dSeqMH;
			m_fMinMass = (float)m_dSeqMH;
			m_State.m_lEqualsS = 0;
			check_parents();
			return true;
		}
// deal with modification at C-terminus, if possible
		if(m_Term.m_bC)	{
			m_Term.m_lState = 2;
			m_dSeqMH += m_pSeqUtilFrag->m_pdAaMod[']'];
			m_Term.m_lN = 0;
			m_Term.m_lC = 1;
			m_State.initialize(m_pSeq,m_lSeqLength);
			m_State.m_dSeqMHS = m_dSeqMH;
			m_fMinMass = (float)m_dSeqMH;
			m_State.m_lEqualsS = 0;
			check_parents();
			return true;
		}
		return false;
	}
// deal with modification at C-terminus, if the N-terminus is in the modified state
	else if(m_Term.m_lState == 1)	{
		if(!m_Term.m_bC)	{
			m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMod['['];
			m_Term.initialize(m_seqUtil.m_pdAaMod['['],m_seqUtil.m_pdAaMod[']']);
			return false;
		}
		m_Term.m_lN = 0;
		m_Term.m_lC = 1;
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMod['['];
		m_dSeqMH += m_pSeqUtilFrag->m_pdAaMod[']'];
		m_Term.m_lState = 2;
		m_State.initialize(m_pSeq,m_lSeqLength);
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float) m_dSeqMH;
		m_State.m_lEqualsS = 0;
		check_parents();
		return true;
	}
// deal with modification at both N- and C- terminus, if possible
	else if(m_Term.m_lState == 2)	{
		if(!m_Term.m_bN)	{
			m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMod[']'];
			m_Term.initialize(m_pSeqUtilFrag->m_pdAaMod['['],m_pSeqUtilFrag->m_pdAaMod[']']);
			return false;
		}
		m_dSeqMH += m_pSeqUtilFrag->m_pdAaMod['['];
		m_Term.m_lN = 1;
		m_Term.m_lC = 1;
		m_State.initialize(m_pSeq,m_lSeqLength);
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float)m_dSeqMH;
		m_Term.m_lState = 3;
		m_State.m_lEqualsS = 0;
		check_parents();
	}
	return bReturn;
}
/*
 * load_next_pam alters the current peptide sequence by iterating through all possible
 * single amino acid polymorphisms. Polymorphisms that are difficult to distinguish by
 * mass spectrometry are skipped to avoid reporting spurious assignments.
 */
bool mscore::load_next_pam(void)
{
	bool bReturn = false;
	if(m_Pam.m_tCount != 0)	{
		m_Pam.m_tAa++;
	}
	m_Pam.m_tCount++;
// Shift to next residue if all possibilities have been checked
	if(m_Pam.m_tAa >= m_Pam.m_tAaTotal)	{
		m_pSeq[m_Pam.m_tPos] = m_Pam.m_pSeqTrue[m_Pam.m_tPos];
		m_Pam.m_tPos++;
		m_Pam.m_tAa = 0;
	}
// Skip mutations indistinguishable from the sequence
	while(m_Pam.m_tPos < m_Pam.m_tLength && check_pam_mass())	{
		if(m_Pam.m_tAa == m_Pam.m_tAaTotal - 1)	{
			m_pSeq[m_Pam.m_tPos] = m_Pam.m_pSeqTrue[m_Pam.m_tPos];
			m_Pam.m_tPos++;
			m_Pam.m_tAa = 0;
		}
		else	{
			m_Pam.m_tAa++;
		}
	}
// Return false if all residues have been checked
	if(m_Pam.m_tPos >= m_Pam.m_tLength)	{
		strcpy(m_pSeq,m_Pam.m_pSeqTrue);
		m_dSeqMH = m_Pam.m_fSeqTrue;
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float)m_dSeqMH;
		m_State.m_lEqualsS = 0;
		check_parents();
		m_Pam.m_tCount = 0;
		return false;
	}
	else	{
		strcpy(m_pSeq,m_Pam.m_pSeqTrue);
		m_dSeqMH = m_Pam.m_fSeqTrue;
		m_dSeqMH += m_pSeqUtilFrag->m_pdAaMass[m_Pam.m_pAa[m_Pam.m_tAa]];
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMass[m_Pam.m_pSeqTrue[m_Pam.m_tPos]];
		m_dSeqMH += m_pSeqUtilFrag->m_pdAaFullMod[m_Pam.m_pAa[m_Pam.m_tAa]];
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaFullMod[m_Pam.m_pSeqTrue[m_Pam.m_tPos]];
		m_pSeq[m_Pam.m_tPos] = m_Pam.m_pAa[m_Pam.m_tAa];
		m_State.initialize(m_pSeq,m_lSeqLength);
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float)m_dSeqMH;
		m_State.m_lEqualsS = 0;
		check_parents();
		return true;
	}
	return bReturn;
}
/*
 * this method checks for modifications that result from trivial changes in mass assignment, e.g.
 * M->F (+16), even though M+16 is being checked a potential modification. getting this function
 * right is critical to the point mutation assignment working properly.
 */
__inline__ bool mscore::check_pam_mass()
{
	const char cTrue = m_Pam.m_pSeqTrue[m_Pam.m_tPos];
	const char cNew = m_Pam.m_pAa[m_Pam.m_tAa];
	const float fTrue = m_pSeqUtilFrag->m_pfAaMass[cTrue] + (float)m_pSeqUtilFrag->m_pdAaFullMod[cTrue];
	const float fNew = m_pSeqUtilFrag->m_pfAaMass[cNew] + (float)m_pSeqUtilFrag->m_pdAaFullMod[cNew];
	if(fabs(fTrue - fNew) < m_fHomoError)	{
		return true;
	}
	if(fabs(fTrue + m_pSeqUtilFrag->m_pdAaMod[cTrue+32] - fNew) < m_fHomoError)	{
		return true;
	}
	if(fabs(fNew + m_pSeqUtilFrag->m_pdAaMod[cNew+32] - fTrue) < m_fHomoError)	{
		return true;
	}
	if(fabs(fTrue + m_pSeqUtilFrag->m_pdAaMod[cTrue+32] - fNew) < m_fHomoError)	{
		return true;
	}
	if(fabs(fNew + m_pSeqUtilFrag->m_pdAaMod[cNew+32] - fTrue - m_pSeqUtilFrag->m_pdAaMod[cTrue+32]) < m_fHomoError)	{
		return true;
	}
	return false;
}

/*
 * load_next_sap alters the current peptide sequence by adding single amino acid
 * polymorphisms that have been annotated in the SAP files specified in the taxonomy
 * file settings.
 */
bool mscore::load_next_sap(void)
{
	bool bReturn = false;
	if(!m_Sap.m_bOk)	{
		return bReturn;
	}
// Return false if all residues have been checked
	if(!m_Sap.next())	{
		memcpy(m_pSeq,m_Sap.m_pSeqTrue,m_lSeqLength);
		m_dSeqMH = m_Sap.m_fSeqTrue;
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float)m_dSeqMH;
		m_State.m_lEqualsS = 0;
		check_parents();
		m_Sap.m_tCount = 0;
		return false;
	}
	else	{
		memcpy(m_pSeq,m_Sap.m_pSeqTrue,m_lSeqLength);
		m_dSeqMH = m_Sap.m_fSeqTrue;
		m_dSeqMH += m_pSeqUtilFrag->m_pdAaMass[m_Sap.m_cCurrent];
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaMass[m_Sap.m_pSeqTrue[m_Sap.m_tPos]];
		m_dSeqMH += m_pSeqUtilFrag->m_pdAaFullMod[m_Sap.m_cCurrent];
		m_dSeqMH -= m_pSeqUtilFrag->m_pdAaFullMod[m_Sap.m_pSeqTrue[m_Sap.m_tPos]];
		m_pSeq[m_Sap.m_tPos] = m_Sap.m_cCurrent;
		m_State.initialize(m_pSeq,m_lSeqLength);
		m_State.m_dSeqMHS = m_dSeqMH;
		m_fMinMass = (float)m_dSeqMH;
		m_State.m_lEqualsS = 0;
		check_parents();
		return true;
	}
	return bReturn;
}
/*
 * load_state is used access the potential sequence modification state machine.
 * it runs the state machine until either:
 *   - a modified sequence that has a mass that is within error of 
 *     one of the mspectrum objects has been found, or
 *   - the last of the modified sequences has been returned
 * it returns true if there are more sequences to consider.
 * NOTE: the last sequence generated by the state machine is the one without any
 *       of the potential modifications in place. maintaining this order is important
 *       to the correct functioning of the logic using in mprocess.
 */
bool mscore::load_state(void)
{
	bool bReturn = run_state_machine();
	if(m_dSeqMH < m_fMinMass)	{
		m_fMinMass = (float)m_dSeqMH;
	}
	while(bReturn)	{
		if(m_State.m_bIsPossible && check_parents())	{
			return bReturn;
		}
		bReturn = run_state_machine();
		if(m_dSeqMH < m_fMinMass)	{
			m_fMinMass =(float) m_dSeqMH;
		}
	}
	return bReturn;
}
/*
 * load_seq is used to call the appropriate add_? method to create the
 * ion type specific scoring array. These arrays are padded to mod(4), so
 * that an optimizing compiler can utilize MMX arrays in the dot method.
 */
bool mscore::load_seq(const unsigned long _t,const long _c)
{
	bool bReturn = true;

	if(T_Y & _t)	{
		bReturn = add_Y(_t,_c);
	}
	else if(T_X & _t)	{
		bReturn = add_X(_t,_c);
	}
	else if(T_A & _t)	{
		bReturn = add_A(_t,_c);
	}
	else if(T_B & _t)	{
		bReturn = add_B(_t,_c);
	}
	else if(T_C & _t)	{
		bReturn = add_C(_t,_c);
	}
	else if(T_Z & _t)	{
		bReturn = add_Z(_t,_c);
	}
	return bReturn;
}

/*
 * modifications have exponential impact on number of sequences that require scoring,
 * so a limit is imposed to keep things like a modification on M, and a sequence with
 * 30 M's from costing too much.
 */

unsigned long mscorestate::M_lMaxModStates = 1 << 12;

/*
 * the state machine for determining all potentially modified residues in a peptide sequence
 * and then generating and scoring them efficiently relies upon run_state_machine.
 *
 * this version calculates the various combinations in order of the number of actual
 * modifications present.  for example, if there are 10 possible modification sites,
 * then all combinations with 1 modification are checked first, then all combinations
 * with 2 modifications, with 3, with 4, etc.
 *
 * this makes it possible to impose a maximum number of states checked, while still
 * checking the most likely combinations first.  i.e. checking all possible combinations
 * of 1-10 modifications on 30 sites is probably better than all combinations of the
 * first 10 sites, as with FIX1.
 */

bool mscore::run_state_machine(void)
{
/*
 * return false if there are no more states
 */
   	 m_State.m_bIsPossible = true;
    if(!m_State.m_bStateS)
      {
			memcpy(m_pSeq,m_State.m_pSeqS,m_lSeqLength);
			m_dSeqMH = m_State.m_dSeqMHS;
            return false;
      }
      else if (m_State.m_lStates >= mscorestate::M_lMaxModStates)
      {
			memcpy(m_pSeq,m_State.m_pSeqS,m_lSeqLength);
			m_dSeqMH = m_State.m_dSeqMHS;
            m_State.m_bStateS = false;
			m_State.m_lStates++;
			return true;
      }
      m_State.m_lStates++;
      bool bReturn = m_State.m_bStateS;
/*
 * deal efficiently with protein modeling sessions that do not require potential modifications
 */
      if(!m_seqUtil.m_bPotential)   {
            m_State.m_bStateS = false;
            return bReturn;
      }
      if(m_State.m_lLastS == 0)     {
            m_State.m_bStateS = false;
            return bReturn;
      }
/*
 * adjust the states of the modification positions
 */
      if(m_State.m_lFilledS > 0 &&
                  m_State.m_piMods[m_State.m_lCursorS] < m_State.m_lLastS - m_State.m_lFilledS + m_State.m_lCursorS)  {
            /*
             * shift current mod through all valid positions.
             */
            m_State.m_piMods[m_State.m_lCursorS]++;
      }
      else if (m_State.m_lCursorS > 0) {
            /*
             * shift the mod to the left of the current mod.
             */
            m_State.m_lCursorS--;
            m_State.m_piMods[m_State.m_lCursorS]++;
            /*
             * if there it is not at its largest possible position, shift all mods to the
             * right back, and start over with them.
             */
            if (m_State.m_piMods[m_State.m_lCursorS] < m_State.m_lLastS - m_State.m_lFilledS + m_State.m_lCursorS)
            {
                  for (unsigned long i = 1; i < m_State.m_lFilledS - m_State.m_lCursorS; i++)
                        m_State.m_piMods[m_State.m_lCursorS + i] = m_State.m_piMods[m_State.m_lCursorS] + i;
                  m_State.m_lCursorS = m_State.m_lFilledS - 1;
            }
      }
      else if (m_State.m_lFilledS < m_State.m_lLastS) {
            /*
             * introduce more possible modifications, and start over calculating all
             * possible combinations.
             */
           m_State.m_lFilledS++;
            if (m_State.m_lFilledS < m_State.m_lLastS)
                  m_State.m_lCursorS = m_State.m_lFilledS - 1;
            for (unsigned long i = 0; i < m_State.m_lFilledS; i++)
                  m_State.m_piMods[i] = i;
	  }
      else {
            /*
             * last state is no modifications.
             */
            m_State.m_lFilledS = 0;
      }
/*
 * reset the sequence string and initial mass
 */
	  memcpy(m_pSeq,m_State.m_pSeqS,m_lSeqLength);
      m_dSeqMH = m_State.m_dSeqMHS;
/*
 * unmodified is the last state
 */
      if (m_State.m_lFilledS == 0)
            m_State.m_bStateS = false;
/*
 * otherwise, use the modification state indices to correctly set the modifications
 * in the sequence string, and add to the parent mass.
 */
      else
      {
/*            for (unsigned long i = 0; i < m_State.m_lFilledS; i++)
            {
                  unsigned long pos = m_State.m_piMods[i];
                  *(m_State.m_ppModsS[pos]) += 32;
                  m_dSeqMH += m_seqUtil.m_pdAaMod[*(m_State.m_ppModsS[pos])]; 
            } */
			m_plAA['s'] = 0;
			m_plAA['t'] = 0;
			m_plAA['n'] = 0;
			m_plAA['q'] = 0;
			m_plAA['y'] = 0;
			unsigned long pos = 0;
			size_t tPos = 0;
			double *pAaMod = m_seqUtil.m_pdAaMod;
			for (unsigned long i = 0; i < m_State.m_lFilledS; i++)	{
				pos = m_State.m_piMods[i];
				*(m_State.m_ppModsS[pos]) += 32;
				tPos = (size_t) *(m_State.m_ppModsS[pos]);
				m_plAA[tPos]++;
				m_dSeqMH += pAaMod[tPos];
			}
			m_State.m_bIsPossible = !(m_plAA['s'] + m_plAA['t'] + m_plAA['y'] > 3 || m_plAA['n'] + m_plAA['q'] > 3);

      }
     return true;
}

/*
 * the score method is called externally to score a loaded peptide against one of the
 * loaded mass spectra. the mass spectrum is refered to by its index number in the
 * m_vSpec mspectrum vector. the sequence has already been loaded via set_seq and/or
 * add_seq.
 */
float mscore::score(const size_t _i)
{
	m_fScore = -1.0;
	m_fHyper = -1.0;
	double dFactor = 1.0;
/*
 * return -1000.0 if there is no sequence available
 */
	if(m_pSeq == NULL)
		return -1000.0;
/*
 * initialize values for the protein modeling session
 */
	prescore(_i);

	double dScore = (float)0.0;
	double dValue = (float)0.0;

	unsigned int lType = T_Y;
	unsigned long lValue = 0;
	unsigned long lValueTotal = 0;
	unsigned long lS = S_Y;

	long lChargeLimit = (long)m_vSpec[m_lId].m_fZ;
	unsigned int uiType = m_vSpec[m_lId].m_uiType;
	//if uiType has not been set, then use m_lType value
	if(uiType == 0)	{
		uiType = (unsigned int)m_lType;
		m_vSpec[m_lId].m_uiType = uiType;
	}
	// 2006.12.01: to reduce problems with high charge states, the hyperscore is
	// calculated on the 2 best fragment ion charge states. in the
	// previous versions the hyperscore used all available charge states, giving
	// a considerable advantage to highly charged parent ions
	vector<double> vFactor;
	long a = 0;
	while(a < lChargeLimit+1)	{
		vFactor.push_back(1.0);
		a++;
	}
	if(lChargeLimit == 1)	{
		lChargeLimit = 2;
	}
	if((uiType & T_C) || (uiType & T_Z))	{
		if(lChargeLimit > 2)	{
			lChargeLimit--;
		}
	}
/*
 * iterate through all of the possible values of the mscore_type enum
 * comparing them against uiType
 */
	while(lType < uiType+1)	{
		lValueTotal = 0;
		dValue = 0.0;
		if(lType & uiType)	{
			a = 1;
			while(a < lChargeLimit)	{
/*
 * load the sequence arrays for each possible charge states for the selected spectrum
 */
				load_seq(lType,a);
				lValue = 0;
/*
 * perform a dot product on each charge state
 */
				dValue += dot(&lValue);
				if(a == 1 && (T_Y & lType) && (long)m_vSpec[m_lId].m_fZ == 2)	{
					unsigned long lTemp = 0;
					add_Y(0,2);
					dValue += dot(&lTemp);
					lValue += lTemp;
				}
				lValueTotal += lValue;
				vFactor[a] *= hfactor(lValue);
				a++;
			}
			dScore += dValue;
		}
		m_pfScore[lS] = (float) dValue;
		m_plCount[lS] = lValueTotal;
/*
 * move on to next value in the mstate_type enum
 */
		lS++;
		lType *= 2;
	}
	dScore *= sfactor();
	m_fScore = (float)dScore;
	// only use the 2 best component hyperscores
	sort(vFactor.begin(),vFactor.end());
	a = (long)vFactor.size()-1;
	if(a > 1)	{
		dFactor = vFactor[a] * vFactor[a-1];
	}
	else	{
		dFactor = vFactor[a];
	}
	dFactor *= dScore;
	if(dFactor > FLT_MAX)	{
		m_fHyper = FLT_MAX;
	}
	else	{
		m_fHyper = (float)dFactor;
	}
/*
 * returning 1.0 for a zero score makes the logic in mprocess easier. see mprocess:create_score
 * to see why.
 */
	if(dScore == 0.0)	{
		dScore = 1.0;
	}
	char *pS = strchr(m_pSeq,'s');
	char *pT = strchr(m_pSeq,'t');

	if(m_bPhosphoBias && m_fHyper < FLT_MAX && (pS || pT))	{
		int iST = 0;
		char *pV = strstr(m_pSeq,"sP");
		while(pV)	{
			iST++;
			pV++;
			pV = strstr(pV,"sP");
		}
		pV = strstr(m_pSeq,"tP");
		while(pV)	{
			iST++;
			pV++;
			pV = strstr(pV,"tP");
		}
		double dV = (double)m_fHyper*(1.0 + 0.001*iST);
		if(dV < FLT_MAX)	{
			 m_fHyper = (float)dV;
		}
		double dNeutral = 0.0;
		unsigned long lNeutral = 0;
		m_dWE = (double)(m_fWidth/m_fErr);
		if((pS && m_pSeqUtilFrag->m_bPhosphoSerine) || (pT && m_pSeqUtilFrag->m_bPhosphoThreonine))	{
			nMap::iterator itMap = m_seqUtil.m_mapNeutralLoss.find(m_lId);
			if(itMap != m_seqUtil.m_mapNeutralLoss.end())	{
				double dHyper = (double)m_fHyper * (double)itMap->second;
				if(dHyper < FLT_MAX)	{
					m_fHyper = (float)dHyper;
				}
			}
			else	{
				dNeutral = m_dSeqMH - (79.966331 + m_seqUtil.m_dWater) - m_seqUtil.m_dProton;
#ifdef PLUGGABLE_SCORING
				lNeutral = mconvert(dNeutral,(long)m_vSpec[m_lId].m_fZ);
#else
				lNeutral = mconvert(dNeutral,(double)m_vSpec[m_lId].m_fZ);
#endif
				float fV = ion_check(lNeutral,m_lId);
				if(fV >= 20.0)	{
					double dHyper = (double)m_fHyper * (double)1.001;
					if(dHyper < FLT_MAX)	{
						m_fHyper = (float)dHyper;
					}
					m_seqUtil.m_mapNeutralLoss.insert(nMap::value_type(m_lId,(float)1.001));
				}
				else	{
					m_seqUtil.m_mapNeutralLoss.insert(nMap::value_type(m_lId,(float)1.0));
				}
			}
		}

	}
	return (float) dScore;
}
/*
 * get the current value of the m_fSeqMH member variable
 */
double mscore::seq_mh()
{
	return m_dSeqMH;
}
/*
 * set true if sequence is to be checked for all possible point mutations
 */
bool mscore::set_pam(const bool _b)
{
	m_bUsePam = _b;
	m_Pam.m_tCount = 0;
	return m_bUsePam;
}
/*
 * set true if corrections to phosphorylation assignments are to be made
 */
bool mscore::set_phospho_bias(const bool _b)
{
	m_bPhosphoBias = _b;
	return m_bPhosphoBias;
}

/*
 * set true if sequence is to be checked for known single amino acid polymorphisms
 */
bool mscore::set_saps(const bool _b,string &_s)
{
	m_bUseSaps = _b;
	m_Sap.reset_value(_s,_b);
	return m_bUseSaps;
}

bool mscore::set_allowed_saps(string &_s)
{
	m_Sap.allowed(_s);
	return m_bUseSaps;
}
/*
 * set the absolute, zero-based start position of the peptide in the full
 * protein sequence
 */
bool mscore::set_pos(const size_t _t)
{
	m_tSeqPos = _t;
	return true;
}
/*
 * set a new peptide sequence for consideration
 */
unsigned long mscore::set_seq(const char *_s,const bool _n,const bool _c,const unsigned long _l,const int _f)
{
	m_tSeqPos = _f;
	if(_s == NULL)
		return 0;
/*
 * adjust arrays if the sequence is longer that m_lSize
 */
	m_lSeqLength = _l;
	if(m_lSeqLength >= m_lSize-1)	{
		if(m_pfSeq != NULL)
			delete m_pfSeq;
		if(m_plSeq != NULL)
			delete m_plSeq;
		if(m_pSeq != NULL)
			delete m_pSeq;
		m_lSize = m_lSeqLength + 16;
		m_pfSeq = new float[m_lSize];
		m_pSeq = new char[m_lSize];
		m_plSeq = new unsigned long[m_lSize];
	}
/*
 * make a copy of the sequence
 */
	strcpy(m_pSeq,_s);
	m_dSeqMH = 0.0;
	unsigned long a = 0;
	m_bIsC = _c;
	m_bIsN = _n;
	m_State.initialize(m_pSeq,m_lSize);
	m_Term.initialize(m_seqUtil.m_pdAaMod['['],m_seqUtil.m_pdAaMod[']']);
	m_State.m_lLastS = 0;
	map<size_t,size_t>::iterator itValue;
	map<size_t,size_t>::iterator itEnd = m_seqUtil.m_mapMotifs.end();
	SMap::iterator itSeq;
	SMap::iterator itSeqEnd = m_seqUtil.m_mapMods.end();
	a = 0;
	if(m_seqUtil.m_bPotentialMotif)	{
		m_seqUtil.clear_motifs(true);
	}
/*
 * calculate the M+H of the sequence
 */
	m_iCharge = 1;
	while(a < m_lSeqLength)	{
		m_dSeqMH += m_seqUtil.m_pdAaMass[m_pSeq[a]] + m_seqUtil.m_pdAaMod[m_pSeq[a]] + m_seqUtil.m_pdAaFullMod[m_pSeq[a]];
		if(m_seqUtil.m_bSequenceMods)	{
			itSeq = m_seqUtil.m_mapMods.find((unsigned long)(m_tSeqPos+a));
			if(itSeq != itSeqEnd)	{
				m_dSeqMH += itSeq->second;
			}
		}
		if(m_seqUtil.m_pdAaMod[m_pSeq[a]+32] != 0.0)
			m_State.add_mod(m_pSeq+a);
		itValue = m_seqUtil.m_mapMotifs.find(m_tSeqPos+a);
		if(m_seqUtil.m_bPotentialMotif)	{
			itValue = m_seqUtil.m_mapMotifs.find(m_tSeqPos+a);
			if(itValue != itEnd){
				m_State.add_mod(m_pSeq+a);
				m_seqUtil.add_mod(m_pSeq[a],itValue->second);
			}
		}
		a++;
	}
	if(m_seqUtil.m_bPotentialMotif)	{
		m_seqUtil.set_motifs();
	}
	m_dSeqMH += m_seqUtil.m_dProton + m_seqUtil.m_dCleaveN + m_seqUtil.m_dCleaveC;
	if(m_Term.m_lN)	{
		m_dSeqMH += m_seqUtil.m_pdAaMod['['];
	}
	if(m_Term.m_lC)	{
		m_dSeqMH += m_seqUtil.m_pdAaMod[']'];
	}
/*
 * deal with protein terminii
 */
	if(m_bIsC)
		m_dSeqMH += m_seqUtil.m_fCT;
	if(m_bIsN)
		m_dSeqMH += m_seqUtil.m_fNT;
	m_dSeqMH += m_seqUtil.m_pdAaFullMod['['];
	m_dSeqMH += m_seqUtil.m_pdAaFullMod[']'];

/*
 * record the M+H value in the mstate object
 */
	m_State.m_dSeqMHS = m_dSeqMH;
	m_fMinMass = (float)m_dSeqMH;
	if(m_bUsePam)	{
		m_Pam.initialize(m_pSeq,m_lSize,(float)m_dSeqMH);
	}
	if(m_bUseSaps)	{
		m_Sap.initialize(m_pSeq,(size_t)m_lSize,(float)m_dSeqMH,_f);
	}
	return m_lSeqLength;
}
/*
 * sets the types of ions (Biemann notation) to be used in scoring a peptide sequence
 */
unsigned long mscore::set_type(const unsigned long _t)
{
	m_lType = _t;
	return m_lType;
}
/*
 * sets the interpretation of the parent ion and fragment ion errors to be either absolute
 * in Daltons, or relative, in parts-per-million. the allowed values are constructed
 * from the values in the enum mscore_error.
 */
unsigned long mscore::set_error(const unsigned long _t)
{
	m_lErrorType = _t;
	return m_lErrorType;
}
/*
 * sets the interpretation of the parent ion and fragment ion errors to be either absolute
 * in Daltons, or relative, in parts-per-million. the allowed values are constructed
 * from the values in the enum mscore_error.
 */
float mscore::set_homo_error(const float _f)
{
	m_fHomoError = _f;
	return m_fHomoError;
}
/*
 * sets the value for the fragment error member value m_fErr. this value is interpreted in two ways:
 * CASE 1: m_lErrorType & T_FRAGMENT_DALTONS is true
 *         in this case, the error is the error in Daltons
 * CASE 2: m_lErrorType & T_FRAGMENT_PPM is true
 *         in this case, the error in Daltons is calculated at m/z = 200.0
 *         this value is used in the blur() method and the width of the
 *         blurred distribution is scaled from the value at m/z = 200.0
 */
float mscore::set_fragment_error(const float _f)
{
	if(_f <= 0.0)
		return 0.0;
	m_fErr = _f;
	if(m_lErrorType & T_FRAGMENT_PPM)	{
		m_fErr = (float)(200.0*m_fErr/1e6);
	}
/*
 * NOTE: the m_fErr value used in the ppm case is: 200 x (the error in ppm)/1000000
 */
	return m_fErr;
}
/*
 * sets the value for the parent error member values m_fParentErrPlus and m_fParentErrMinus.
 */
float mscore::set_parent_error(const float _f,const bool _b)
{
	if(_b)	{
		if(_f < 0.0)
			m_fParentErrPlus = 0.0;
		else
			m_fParentErrPlus = _f;
	}
	else	{
		if(_f < 0.0)
			m_fParentErrPlus = 0.0;
		else
			m_fParentErrMinus = _f;
	}
	return _f;
}
/*
 * sets the value for the m_bIsotopeError member, which checks for isotope peak assignment errors
 */
bool mscore::set_isotope_error(const bool _b)
{
	m_bIsotopeError = _b;
	return m_bIsotopeError;
}

/*
 * sets the mass type for fragment ions
 */
void mscore::set_fragment_masstype(masscalc::massType _t)
{
	if (_t == masscalc::average)	{
		m_seqUtilAvg.set_modified(true);
		m_pSeqUtilFrag = &m_seqUtilAvg;
	}
	else	{
		m_pSeqUtilFrag = &m_seqUtil;
	}
}

/*
 * return true if any parent ion is within error of the current peptide's modified M+H
 * and return the number of spectra that have M+H values greater than or within error
 * of the current peptide's modified M+H
 */
bool mscore::test_parents(size_t &_t)	{
	size_t a = 0;
	const float fSeqMH = (float)m_dSeqMH;
	if(m_lSpectra > 100)	{
		size_t h = m_lSpectra/10;
		while(!(m_vDetails[a] < fSeqMH) && a < 9*h)	{
			a += h;
		}
	}
	const size_t tLimit = (size_t)m_lSpectra;
	while(a < tLimit)	{
		if(m_vDetails[a] == fSeqMH)	{
			_t = tLimit - a;
			return true;
		}
		a++;
	}
	return false;
}

bool mscore::reset_permute()
{
	m_psPermute.m_tPos = 0;
	m_psPermute.m_tEnd = m_lSeqLength-2;
	if(m_lSeqLength > m_psPermute.m_lSize)	{
		delete m_psPermute.m_pPerm;
		delete m_psPermute.m_pSeq;
		m_psPermute.m_lSize = m_lSeqLength + 16;
		m_psPermute.m_pPerm = new char[m_psPermute.m_lSize+1];
		m_psPermute.m_pSeq = new char[m_psPermute.m_lSize+1];
	}
	strcpy(m_psPermute.m_pSeq,m_pSeq);
	m_psPermute.m_bRev = true;
	return true;
}

bool mscore::permute()
{
	if(m_psPermute.m_tPos == m_psPermute.m_tEnd && m_psPermute.m_bRev)	{
		strcpy(m_pSeq,m_psPermute.m_pSeq);
		string strTemp;
		string strSeq = m_pSeq;
		string::reverse_iterator itS = strSeq.rbegin();
		string::reverse_iterator itE = strSeq.rend();
		while(itS != itE)	{
			strTemp += *itS;
			itS++;
		}
		strcpy(m_pSeq,strTemp.c_str());
		m_psPermute.m_bRev = false;
		m_psPermute.m_tPos = 0;
	}
	if(m_psPermute.m_tPos == m_psPermute.m_tEnd)	{
		strcpy(m_pSeq,m_psPermute.m_pSeq);
		return false;
	}
	memcpy(m_psPermute.m_pPerm + 1,m_pSeq,m_lSeqLength);
	m_psPermute.m_pPerm[0] = m_psPermute.m_pPerm[m_lSeqLength];
	m_psPermute.m_pPerm[m_lSeqLength] = '\0';
	memcpy(m_pSeq,m_psPermute.m_pPerm,m_lSeqLength);
	m_psPermute.m_tPos++;
	return true;
}

/*
 * mscoremanager contains static short-cuts for dealing with mscore
 * plug-ins.
 */
const char* mscoremanager::TYPE = "scoring, algorithm";

/*
 * create_mscore creates the correct mscore for the given set of XML
 * parameters.
 */
mscore* mscoremanager::create_mscore(XmlParameter &_x)
{
	string strValue;
	string strKey = TYPE;
	if (!_x.get(strKey,strValue))
		strValue = "tandem";

	return (mscore*) mpluginmanager::get().create_plugin(TYPE, strValue.data());
}

/*
 * register_factory registers a factory with the mpluginmanager for
 * creating mscore derived objects.
 */
void mscoremanager::register_factory(const char* _spec, mpluginfactory* _f)
{
	mpluginmanager::get().register_factory(TYPE, _spec, _f);
}

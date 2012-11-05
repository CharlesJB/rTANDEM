/*
 Copyright (C) 2003 Ronald C Beavis, all rights reserved
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

// File version: 2003-07-01
// File version: 2005-01-01

/*
 * the mspectrumcondition object analyzes tandem mass spectra and alters them to
 * be suitable for searching. it loads a set of parameters from an XmlParameters
 * object and uses those values for the spectrum analysis. this object contains
 * more "get" and "set" type methods than most of the other objects in this
 * project: there was considerable experimentation necessary to refine the spectrum
 * analysis methods and these additional methods were used for these experiments.
  *
 * mspectrumcondition performs a series of operations on a spectrum list to make
 * compatible with the scoring methods. the spectrum can be normalized, adjacent
 * peaks removed, parent ion masses excluded, low mass parents excluded, highly
 * charged parents excluded, low mass fragments excluded, and low intensity
 * fragment ions excluded.
*/
#include "stdafx.h"
#include <algorithm>
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mspectrumcondition.h"
#include "mscore.h"

bool lessThanMImz(const mi &_l,const mi &_r);

/*
 * global less than operator for mi classes: to be used in sort operations to achieve
 * list ordered from lowest to highest m/z
 */
bool lessThanMImz(const mi &_l,const mi &_r)
{
	return _l.m_fM < _r.m_fM;
}

/*
 * mspectrumcondition performs a series of operations on a spectrum list to make
 * compatible with the scoring methods. the spectrum can be normalized, adjacent
 * peaks removed, parent ion masses excluded, low mass parents excluded, highly
 * charged parents excluded, low mass fragments excluded, and low intensity
 * fragment ions excluded.
 */
mspectrumcondition::mspectrumcondition(void)
{
	m_bUseMaxPeaks = true;
	m_bUseDynamicRange = true;
	m_bUseMinMass = true;
	m_bUseParent = true;
	m_bUseLowestMass = true;
	m_bUseMinSize = true;
	m_bUseNoiseSuppression = true;
	m_bUseChargeSuppression = true;
	m_bCondition = true;
	m_bUseNeutralLoss = false;
	m_bUsePhosphoDetection = false;

	m_tMaxPeaks = 50;
	m_fDynamicRange = 100.0;
	m_fMinMass = 500.0;
	m_fLowestMass = 150.0;
	m_fParentLower = 2.0;
	m_fParentUpper = 2.0;
	m_lMinSize = 5;
	m_lMaxCharge = 3;
	m_fNeutralLoss = 0.0;
	m_fNeutralLossWidth = 0.0;
	m_fFactor = 1.0;
	m_fMaxZ = 4.0;
}

mspectrumcondition::~mspectrumcondition(void)
{
}
/*
 * condition is the method called by mprocess to condition a particular mspectrum object.
 * by checking the appropriate flags, the various properties of a spectrum can be altered,
 * and the spectrum rejected if it does not pass the appropriate tests. 
 * NOTE: the order of operations is important here: do not change it unless you feel you
 *       fully understand why things were done in this order.
 */

bool mspectrumcondition::condition(mspectrum &_s, mscore &_score)
{
/*
 * bail out if conditioning has been turned off
 */
	if(m_bUsePhosphoDetection)	{
		if(!find_loss(_s,98.0,3.0))	{
			return false;
		}
	}
	if(m_bUseAllowedNeutralLosses)	{
		size_t m = 0;
		size_t tSize = m_vdNeutralLosses.size();
		bool bOK = false;
		while(m < tSize && !bOK)	{
			bOK = find_loss(_s,(float)m_vdNeutralLosses[m],(float)0.5,(float)0.05);
			m++;
		}
		if(!bOK)	{
			return false;
		}
	}
	_s.m_vMINeutral.clear();
	sort(_s.m_vMI.begin(),_s.m_vMI.end(),lessThanMImz);
	if(_s.m_vdStats.size() == 0)	{
		vector<mi>::iterator itMI = _s.m_vMI.begin();
		vector<mi>::iterator itEnd = _s.m_vMI.end();
		double dSum = 0.0;
		double dMax = 0.0;
		while(itMI != itEnd)	{
			if(dMax < itMI->m_fI)	{
				dMax = itMI->m_fI;
			}
			dSum += itMI->m_fI;
			itMI++;
		}
		_s.m_vdStats.push_back(dSum);
		_s.m_vdStats.push_back(dMax);
		_s.m_vdStats.push_back(m_fFactor);
	}
	if(_s.m_fZ > m_fMaxZ)	{
		return false;
	}
	if(!m_bCondition)	{
		return true;
	}
/*
 * give the score object a chance to perform score specific modifications
 * or veto this spectrum.
 */
	if(!_score.precondition(_s)) {
		return false;
	}
/*
 * check the parent ion charge and return false if it is too large
 */
	if(m_bUseNoiseSuppression)	{
		/*
		 * check the parent ion mass against the minimum mass allowed (about 600.0 is usually safe)
		*/
		if(m_bUseMinMass)	{
			if(_s.m_dMH < m_fMinMass)
				return false;
		}
		if(m_bUseChargeSuppression)	{
			if((long)(_s.m_fZ+0.5) > m_lMaxCharge)	{
				return false;
			}
		}
	}
	size_t tSize = _s.m_vMI.size();
	size_t a = 0;
	float fMaxI = 0;
/*
 * this method doesn't really remove isotopes: it cleans up multiple intensities within one Dalton
 * of each other.
 */
//	sort(_s.m_vMI.begin(),_s.m_vMI.end(),lessThanMImz);
	remove_isotopes(_s);
/*
 * remove ions near the parent ion m/z
 */
	if(m_bUseParent)	{
		remove_parent(_s);
	}
/*
 * remove low mass immonium ions prior to normalization
 */
	if(m_bUseLowestMass)	{
		remove_low_masses(_s);
	}
/*
 * normalize the spectrum 
 */
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	vector<mi>::iterator itEnd = _s.m_vMI.end();
	double dSum = 0.0;
	double dMax = 0.0;
	while(itMI != itEnd)	{
		if(dMax < itMI->m_fI)	{
			dMax = itMI->m_fI;
		}
		dSum += itMI->m_fI;
		itMI++;
	}
	_s.m_vdStats.clear();
	_s.m_vdStats.push_back(dSum);
	_s.m_vdStats.push_back(dMax);
	_s.m_vdStats.push_back(m_fFactor);
	if(m_bUseDynamicRange)	{
		dynamic_range(_s);
	}
	if(m_bUseNeutralLoss)	{
		remove_neutral(_s);
	}
/*
* reject the spectrum if there aren't enough peaks 
*/
	if(m_bUseMinSize)	{
		if((long)_s.m_vMI.size() < m_lMinSize)
			return false;
	}
/*
 * check to see if the spectrum has the characteristics of noise 
 */
	if(m_bUseNoiseSuppression)	{
		if(is_noise(_s))	{
			return false;
		}
	}
/*
* retrieve the N most intense peaks
*/
	clean_isotopes(_s);
	if(m_bUseMaxPeaks)	{
		sort(_s.m_vMI.begin(),_s.m_vMI.end(),lessThanMI);
		remove_small(_s);
	}
	sort(_s.m_vMI.begin(),_s.m_vMI.end(),lessThanMImz);
	itMI = _s.m_vMI.begin();
	itEnd = _s.m_vMI.end();
	dSum = 0.0;
	dMax = 0.0;
	while(itMI != itEnd)	{
		if(dMax < itMI->m_fI)	{
			dMax = itMI->m_fI;
		}
		dSum += itMI->m_fI;
		itMI++;
	}
	_s.m_vdStats[0] = dSum*m_fFactor;
	_s.m_vdStats[1] = dMax*m_fFactor;
	_s.m_vdStats[2] = m_fFactor;
	return true;
}
/*
 * check_neutral looks for the loss of water or ammonia in a spectrum. this
 * method is not used in this implementation: it is for further experimentation
 */
bool mspectrumcondition::find_loss(mspectrum &_s,const float _d,const float _t,const float _p)
{
	sort(_s.m_vMI.begin(),_s.m_vMI.end(),lessThanMI);
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	vector<mi>::iterator itEnd = _s.m_vMI.end();
	float fMax = itMI->m_fI;
	while(itMI != itEnd)	{
		if(itMI->m_fI > fMax)	{
			fMax = itMI->m_fI;
		}
		itMI++;
	}
	fMax *= _p;
	itMI = _s.m_vMI.begin();
	float fD = (float)(1.007276 + ((_s.m_dMH - 1.007276 - _d)/_s.m_fZ));
	while(itMI != itEnd)	{
		if(fabs(itMI->m_fM - fD) <= _t && itMI->m_fI >= fMax)	{
			return true;
		}
		itMI++;
	}
	return false;
}

/*
 * check_neutral looks for the loss of water or ammonia in a spectrum. this
 * method is not used in this implementation: it is for further experimentation
 */
bool mspectrumcondition::check_neutral(mspectrum &_s)
{
	sort(_s.m_vMI.begin(),_s.m_vMI.end(),lessThanMI);
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	vector<mi>::iterator itEnd = _s.m_vMI.end();
	vector<mi>::iterator itLast = itMI;
	unsigned long a = 0;
	unsigned long lCount = 0;
	float fM = 0.0;
	while(a < 10)	{
		itMI = itLast;
		while(itMI != itEnd && itMI->m_fM < 300.0)	{
			itMI++;
		}
		if(itMI == itEnd)
			break;
		fM = (float)(itMI->m_fM - 18.0);
		itMI++;
		itLast = itMI;
		while(itMI < itEnd)	{
			if(fabs(fM - itMI->m_fM) < 2.5)	{
				lCount++;
				break;
			}
			itMI++;
		}
		a++;
	}
	if(lCount < 1)
		return false;
	return true;
}
/*
 * remove isotopes removes multiple entries within 0.95 Da of each other, retaining
 * the highest value. this is necessary because of the behavior of some peak 
 * finding routines in commercial software
 */

bool mspectrumcondition::remove_isotopes(mspectrum &_s)
{
       if(_s.m_vMI.size() < 2)
                return true;

        vector<mi>::iterator itMI_ecriture = _s.m_vMI.begin();
        vector<mi>::iterator itMI_lecture = _s.m_vMI.begin() + 1;
        vector<mi>::const_iterator itMI_fin = _s.m_vMI.end();
		vector<mi> miTemp;
		float fEcriture = itMI_ecriture->m_fM;
        while( itMI_lecture != itMI_fin )
        {
                if( ( itMI_lecture->m_fM - fEcriture) >= 0.95 || itMI_lecture->m_fM < 200.0)
                {
						miTemp.push_back(*itMI_ecriture);
                        itMI_ecriture = itMI_lecture;
						fEcriture = itMI_ecriture->m_fM;
                }
                else if( itMI_lecture->m_fI > itMI_ecriture->m_fI )
                {
                        *itMI_ecriture = *itMI_lecture;
                }
                itMI_lecture++;
        }
		miTemp.push_back(*itMI_ecriture);
		_s.m_vMI = miTemp;
        return true;
}

/*
 * clean_isotopes removes peaks that are probably C13 isotopes
 */
bool mspectrumcondition::clean_isotopes(mspectrum &_s)
{

       if( _s.m_vMI.size() < 2 )
                return true;

        vector<mi>::iterator itMI_ecriture = _s.m_vMI.begin();
        vector<mi>::iterator itMI_lecture = _s.m_vMI.begin() + 1;
        vector<mi>::const_iterator itMI_fin = _s.m_vMI.end();
 		vector<mi> miTemp;
 		float fEcriture = itMI_ecriture->m_fM;
        while( itMI_lecture != itMI_fin )
        {
                if( ( itMI_lecture->m_fM - fEcriture) >= 1.5  || itMI_lecture->m_fM < 200.0)
                {
						miTemp.push_back(*itMI_ecriture);
                        itMI_ecriture = itMI_lecture;
						fEcriture = itMI_ecriture->m_fM;
                }
                else if( itMI_lecture->m_fI > itMI_ecriture->m_fI )
                {
                        *itMI_ecriture = *itMI_lecture;
                }
                itMI_lecture++;
        }
		miTemp.push_back(*itMI_ecriture);
		_s.m_vMI = miTemp;
        return true;

}
/*
 * is_noise attempts to determine if the spectrum is simply noise. if the spectrum
 * does not have any peaks within a window near the parent ion mass, it is considered
 * noise.
 */

bool mspectrumcondition::is_noise(mspectrum &_s)
{
	if(!m_bUseNoiseSuppression)
		return false;
	size_t a = 0;
	size_t tSize = _s.m_vMI.size();
	float fZ = _s.m_fZ;
	float fMax = (float)(_s.m_dMH/fZ);
	if(fZ == 1)	{
		fMax = (float)(_s.m_dMH - 600.0);
	}
	if(fZ == 2)	{
		fMax = (float)(_s.m_dMH - 600.0);
	}
	while(a < tSize)	{
		if(_s.m_vMI[a].m_fM > fMax)
			return false;
		a++;
	}
	return true;
}
/*
 * load makes internal copies of the relavent parameters found in the XmlParameter object
 */
bool mspectrumcondition::load(XmlParameter &_x)
{
	string strKey;
	string strValue;
	strKey = "spectrum, dynamic range";
	m_bUseDynamicRange = _x.get(strKey,strValue);
	if(m_bUseDynamicRange)	{
		m_fDynamicRange = (float)atof(strValue.c_str());
	}
	strKey = "spectrum, total peaks";
	m_bUseMaxPeaks = _x.get(strKey,strValue);
	if(m_bUseMaxPeaks)	{
		m_tMaxPeaks = atoi(strValue.c_str());
	}
	strKey = "spectrum, minimum peaks";
	m_bUseMinSize = _x.get(strKey,strValue);
	if(m_bUseMinSize)	{
		m_lMinSize = atoi(strValue.c_str());
	}
	strKey = "spectrum, minimum parent m+h";
	m_bUseMinMass = _x.get(strKey,strValue);
	if(m_bUseMinMass)	{
		m_fMinMass = (float)atof(strValue.c_str());
	}
	strKey = "spectrum, minimum fragment mz";
	m_bUseLowestMass = _x.get(strKey,strValue);
	if(m_bUseLowestMass)	{
		m_fLowestMass = (float)atof(strValue.c_str());
	}
	strKey = "spectrum, use conditioning";
	bool bCondition = _x.get(strKey,strValue);
	if(bCondition)	{
		if(strValue == "yes")	{
			m_bCondition = true;
		}
		else if(strValue == "no")	{
			m_bCondition = false;
		}
		else	{
			m_bCondition = true;
		}
	}
	strKey = "spectrum, use noise suppression";
	bool bUseNoiseSuppression = _x.get(strKey,strValue);
	if(bUseNoiseSuppression)	{
		if(strValue == "yes")	{
			m_bUseNoiseSuppression = true;
		}
		else if(strValue == "no")	{
			m_bUseNoiseSuppression = false;
		}
		else	{
			m_bUseNoiseSuppression = true;
		}
	}
	strKey = "spectrum, use neutral loss window";
	bool bUseNeutralLoss = _x.get(strKey,strValue);
	if(bUseNeutralLoss)	{
		if(strValue == "yes")	{
			m_bUseNeutralLoss = true;
		}
		else if(strValue == "no")	{
			m_bUseNeutralLoss = false;
		}
		else	{
			m_bUseNeutralLoss = false;
		}
	}
	if(m_bUseNeutralLoss)	{
		strKey = "spectrum, neutral loss window";
		if(_x.get(strKey,strValue))	{
			m_fNeutralLossWidth = (float)atof(strValue.c_str());
		}
		strKey = "spectrum, neutral loss mass";
		if(_x.get(strKey,strValue))	{
			m_fNeutralLoss = (float)atof(strValue.c_str());
		}
	}
	strKey = "spectrum, allowed neutral losses";
	m_bUseAllowedNeutralLosses = _x.get(strKey,strValue);
	if(m_bUseAllowedNeutralLosses)	{
		m_vdNeutralLosses.clear();
		string strV;
		size_t tLength = strValue.size();
		size_t a = 0;
		while(a < tLength)	{
			if(!isspace(strValue[a]))	{
				strV += strValue[a];
			}
			a++;
		}
		tLength = strV.size();
		if(tLength > 0)	{
			a = 0;
			size_t b = strV.find(',',a);
			double dValue = 0.0;
			while(b != strV.npos)	{
				dValue = atof(strV.substr(a,b-a).c_str());
				if(dValue != 0.0)	{
					m_vdNeutralLosses.push_back(dValue);
				}
				a = b + 1;
				b = strV.find(',',a);
			}
			dValue = atof(strV.substr(a,tLength-a).c_str());
			if(dValue != 0.0)	{
				m_vdNeutralLosses.push_back(dValue);
			}
			if(m_vdNeutralLosses.empty())	{
				m_bUseAllowedNeutralLosses = false;
			}
		}
		else	{
				m_bUseAllowedNeutralLosses = false;
		}
	}
	strKey = "spectrum, maximum parent charge";
	if(_x.get(strKey,strValue))	{
		m_fMaxZ = (float)atof(strValue.c_str());
		if(m_fMaxZ < 1.0)	{
			m_fMaxZ = 4.0;
		}
	}
	return true;
}
/*
 * use the dynamic range parameter to set the maximum intensity value for the spectrum.
 * then remove all peaks with a normalized intensity < 1
 */
bool mspectrumcondition::dynamic_range(mspectrum &_s)
{
	if (!m_bUseDynamicRange) {
		return false;
	}
	size_t tSize = _s.m_vMI.size();
	float fI = 1.0;
	if(tSize > 0)	{
		fI = _s.m_vMI[0].m_fI;
	}
	size_t a = 0;
	while(a < tSize)	{
		if(_s.m_vMI[a].m_fI > fI)	{
			fI = _s.m_vMI[a].m_fI;
		}
		a++;
	}
	m_fFactor = fI/m_fDynamicRange;
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	while(itMI != _s.m_vMI.end())	{
		itMI->m_fI /= m_fFactor;
		if(itMI->m_fI < 1.0)	{
			itMI = _s.m_vMI.erase(itMI);
		}
		else	{
			itMI++;
		}
	}
	return true;
}
/*
 * limit the total number of peaks used
 */
bool mspectrumcondition::remove_small(mspectrum &_s)
{
	if(!m_bUseMaxPeaks)	{
		return false;
	}
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	// The test line was changed in 2006.09.01 to correct a problem
	// created by VC++ 2005. This new version uses a strict bounds checking
	// style for STL iterators, that cause a run time error if an iterator
	// longer than .end() is produced by incrementing the iterator. The previous
	// test was 
	//
	//				if(itMI+m_tMaxPeaks < _s.m_vMI.end())
	//
	// which caused the error when the statement was false.
	if(m_tMaxPeaks < _s.m_vMI.size())	{
		_s.m_vMI.erase(itMI+m_tMaxPeaks,_s.m_vMI.end());
	}
	return true;
}
/*
 * set up m/z regions to ignore: those immediately below the m/z of the parent ion
 * which will contain uninformative neutral loss ions, and those immediately above
 * the parent ion m/z, which will contain the parent ion and its isotope pattern
 */
bool mspectrumcondition::remove_parent(mspectrum &_s)
{
	if(!m_bUseParent)
		return false;
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	float fParentMz = (float)(1.00727+(_s.m_dMH - 1.00727)/_s.m_fZ);
	while(itMI != _s.m_vMI.end())	{
		if(fParentMz - itMI->m_fM >= 0.0 && fParentMz - itMI->m_fM < m_fParentLower/_s.m_fZ)	{
			itMI = _s.m_vMI.erase(itMI);
		}
		else if(itMI->m_fM - fParentMz > 0.0 && itMI->m_fM - fParentMz < m_fParentUpper/_s.m_fZ)	{
			itMI = _s.m_vMI.erase(itMI);
		}
		else	{
			itMI++;
		}
	}
	return true;
}
/*
 * remove_low_masses deletes peaks with m/z values below the m_fLowestMass member value
 */
bool mspectrumcondition::remove_low_masses(mspectrum &_s)
{
	if(!m_bUseLowestMass)
		return false;
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	while(itMI != _s.m_vMI.end())	{
		if(itMI->m_fM > m_fLowestMass)	{
			break;
		}
		itMI++;
	}
	_s.m_vMI.erase(_s.m_vMI.begin(),itMI);
	return true;
}
/*
 * remove_neutral
 */
bool mspectrumcondition::remove_neutral(mspectrum &_s)
{
	if(!m_bUseNeutralLoss)
		return false;
	vector<mi>::iterator itMI = _s.m_vMI.begin();
	while(itMI != _s.m_vMI.end())	{
		if(fabs((_s.m_dMH - itMI->m_fM) - m_fNeutralLoss) <= m_fNeutralLossWidth)	{
			_s.m_vMINeutral.push_back(*itMI);
			itMI = _s.m_vMI.erase(itMI);
		}
		else	{
			itMI++;
		}
	}
	return true;
}
/*
 * get_noise_suppression gets the current value of the m_bUseNoiseSuppression member
 */
bool mspectrumcondition::get_noise_suppression(void)
{
	return m_bUseNoiseSuppression;
}
/*
 * set_max_peaks sets the current value of the m_tMaxPeaks member
 */
bool mspectrumcondition::set_max_peaks(const long _p)
{
	m_tMaxPeaks = _p;
	return true;
}
/*
 * set_min_size sets the current value of the m_lMinSize member
 */
bool mspectrumcondition::set_min_size(const long _p)
{
	m_lMinSize = _p;
	return true;
}
/*
 * set_dynamic_range sets the current value of the m_fDynamicRange member
 */
bool mspectrumcondition::set_dynamic_range(const float _p)
{
	m_fDynamicRange = _p;
	return true;
}
/*
 * set_parent_exclusion sets the current value of the m_fParentLower & m_fParentUpper members
 */
bool mspectrumcondition::set_parent_exclusion(const float _l,const float _u)
{
	m_fParentLower = _l;
	m_fParentUpper = _u;
	return true;
}
/*
 * set_min_mass sets the current value of the m_fMinMass member
 */
bool mspectrumcondition::set_min_mass(const float _m)
{
	m_fMinMass = _m;
	return true;
}
/*
 * set_max_charge sets the current value of the m_lMaxCharge member
 */
bool mspectrumcondition::set_max_charge(const long _m)
{
	m_lMaxCharge = _m;
	return true;
}
/*
 * set_lowest_mass sets the current value of the m_fLowestMass member
 */
bool mspectrumcondition::set_lowest_mass(const float _m)
{
	m_fLowestMass = _m;
	return true;
}
/*
 * use_condition sets the current value of the m_bCondition member
 */
bool mspectrumcondition::use_condition(const bool _f)
{
	m_bCondition = _f;
	return m_bCondition;
}
/*
 * use_min_size sets the current value of the m_bUseMinSize member
 */
bool mspectrumcondition::use_min_size(const bool _f)
{
	m_bUseMinSize = _f;
	return m_bUseMinSize;
}
/*
 * use_max_peaks sets the current value of the m_bUseMaxPeaks member
 */
bool mspectrumcondition::use_max_peaks(const bool _f)
{
	m_bUseMaxPeaks = _f;
	return m_bUseMaxPeaks;
}
/*
 * use_charge_suppression sets the current value of the m_bUseChargeSuppression member
 */
bool mspectrumcondition::use_charge_suppression(const bool _f)
{
	m_bUseChargeSuppression = _f;
	return m_bUseChargeSuppression;
}
/*
 * use_dynamic_range sets the current value of the m_bUseDynamicRange member
 */
bool mspectrumcondition::use_dynamic_range(const bool _f)
{
	m_bUseDynamicRange = _f;
	return m_bUseDynamicRange;
}
/*
 * use_noise_suppression sets the current value of the m_bUseNoiseSuppression member
 */
bool mspectrumcondition::use_noise_suppression(const bool _f)
{
	m_bUseNoiseSuppression = _f;
	return m_bUseNoiseSuppression;
}
/*
 * use_lowest_mass sets the current value of the m_bUseLowestMass member
 */
bool mspectrumcondition::use_lowest_mass(const bool _f)
{
	m_bUseLowestMass = _f;
	return m_bUseLowestMass;
}
/*
 * use_min_mass sets the current value of the m_bUseMinMass member
 */
bool mspectrumcondition::use_min_mass(const bool _f)
{
	m_bUseMinMass = _f;
	return m_bUseMinMass;
}
/*
 * use_parent_exclusion sets the current value of the m_bUseParent member
 */
bool mspectrumcondition::use_parent_exclusion(const bool _f)
{
	m_bUseParent = _f;
	return m_bUseParent;
}


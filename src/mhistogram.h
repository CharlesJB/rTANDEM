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

#ifndef MHISTOGRAM_H
#define MHISTOGRAM_H

#include <string.h>  // rTANDEM
// File version: 2003-07-01

/*
 *	The mhistogram class is used for storing information about scoring in the mspectrum class.
 * Two types of histogram classes are used: the mhistogram which stores information about
 * convolution and hyper scores and the count_mhistogram class which is used to store the number
 * of specific ion types (e.g. y or b ions). The mhistgram class converts the scores into base-10 logs
 * to make the histogram more compact and to linearize the extreme value distribution, which
 * governs the scoring distribution. In order to calculate expectation values from these distributions,
 * they are first converted to survival functions (by the survival methdod), the outlying non-random 
 * scores are removed from the distribution and the high scoring tail of the distribution is fitted
 * to a log-log linear distribution using least-squares in the model method.
 * mhistogram is included in the mspectrum.h file only
 */

class mhistogram	{
public:
	mhistogram(void) {
		init();
	}

	mhistogram(const mhistogram& rhs) {
		init();
		(*this) = rhs;
	}
	unsigned long m_ulCount;
	void init() {
		m_ulCount = 0;
		m_lLength = 0;
		m_pList = NULL;
		m_fA0 = (float)4.8;
		m_fA1 = (float)-0.28;
		m_dProteinFactor = 1.0;
		m_lSum = 0;
		m_dLimit = 1.0e-15;
	}

	virtual ~mhistogram(void) {
		delete [] m_pList;
	}
protected:
	double m_dProteinFactor; // a weighting factor
	float m_fA0; // the intercept of the least-squares fit performed in model
	float m_fA1; // the slope of the least-squares fit performed in model
	long m_lLength; // the length of the histogram
private:
	vector<long> m_vlSurvive; // the survival function array, reduced to [96] from [256]
	unsigned short* m_pList; // doubling buffer
	long m_lSum;
	double m_dLimit;

public:

/*
 * expect uses the equation derived in model to convert scores into expectation values
 * survival and model must be called before expect will return reasonable values
 */
	float expect(const float _f)	{
		double dR = pow(10.0,(double)(m_fA0 + m_fA1*_f));
		if(dR < m_dLimit)	{
			dR = m_dLimit;
		}
		return (float)dR;
	}
/*
 * list returns a specified value from the stochastic distribution of the histogram
 */
	long list(const unsigned long _l)	{
		if (_l >= (unsigned long) m_lLength)
			return 0;
		return m_pList[_l];
	}
/*
 * survive returns a specified value from the stochastic distribution of the survival function
 */
	long survive(const unsigned long _l)	{
		return m_vlSurvive[_l];
	}
	long sum()	{
		return m_lSum;
	}
/*
 * length returns the maximum length of the histogram
 */
	long length(void)	{
		return m_lLength;
	}
/*
 * a0 returns the first parameter in the linear fit to the survival function
 */
	float a0(void)	{
		return m_fA0;
	}
/*
 * a1 returns the first parameter in the linear fit to the survival function
 */
	float a1(void)	{
		return m_fA1;
	}
/*
 * expect_protein allows for a modified expectation value to be used for proteins
 */
	float expect_protein(const float _f)	{
		double dR = pow(10.0,(double)(m_fA0 + m_fA1*_f))*m_dProteinFactor;
		if(dR < m_dLimit)	{
			dR = m_dLimit;
		}
		return (float)dR;
	}
/*
 * set_protein_factor simply sets the value for m_dProteinFactor
 */
	bool set_protein_factor(const double _d)	{
		if(_d <= 0.0)
			return false;
		m_dProteinFactor = _d;
		return true;
	}
/*
 * reset zeros the histogram array m_pList
 */
	virtual bool clear()	{
		long a = 0;
		while(a < m_lLength)	{
			m_pList[a] = 0;
			a++;
		}
		m_ulCount = 0;
		return true;
	}

/*
 * add increments the appropriate value in m_pList for a score value
 */
	virtual long add(const float _f)	{
		long lValue = (long)(_f + 0.5);
		// If buffer too small, double its size.
		// Buffer must have empty bucket at end to output
		// histogram correctly for reports.
		if(lValue > m_lLength - 2)	{
			long lLengthNew = lValue + 2;

			unsigned short* pListNew = new unsigned short[lLengthNew];
			memset(pListNew, 0, lLengthNew * sizeof(unsigned short));

			if (m_pList != NULL) {
				memcpy(pListNew, m_pList, m_lLength * sizeof(unsigned short));
				delete [] m_pList;
			}

			m_pList = pListNew;
			m_lLength = lLengthNew;
		}
		
		if(m_pList[lValue] < 0xFFFE)	{
			m_pList[lValue]++;
		}
		m_ulCount++;
		return lValue;
	}
/*
 * survival converts the scoring histogram in m_pList into a survival function
 */
	bool survival()
	{
		if (m_lLength == 0)
			return false;

		long a = m_lLength - 1;
		long lSum = 0;
		long *plValues = new long[m_lLength];
/*
 * first calculate the raw survival function
 */
		while(a > -1)	{
			lSum += m_pList[a];
			plValues[a] = lSum;
			a--;
		}
		a = 0;
/*
 * determine the value of the survival function that is 1/5 of the maximum
 */
		const long lPos = plValues[0]/5;
		while(a < m_lLength && plValues[a] > lPos)	{
			a++;
		}
		const long lMid = a;
		a = m_lLength - 1;
		while(a > -1 && plValues[a] == 0)	{
			a--;
		}
		lSum = 0;
		long lValue = 0;
/*
 * remove potentially valid scores from the stochastic distribution
 */
		while(a > 0)	{
			if(plValues[a] == plValues[a-1] 
					&& plValues[a] != plValues[0]
					&& a > lMid)	{
				lSum = plValues[a];
				lValue = plValues[a];
				a--;
				while(plValues[a] == lValue)	{
					plValues[a] -= lSum;
					a--;
				}
			}
			else	{
				plValues[a] -= lSum;
				a--;
			}
		}
		plValues[a] -= lSum;
		a = 0;
/*
 * replace the scoring distribution with the survival function
 */
		m_vlSurvive.clear();
		while(a < m_lLength)	{
			m_vlSurvive.push_back(plValues[a]);
			a++;
		}
		delete plValues;
		m_lSum = m_vlSurvive[0];
		return true;
	}

	bool clear_survive()
	{
		m_vlSurvive.clear();
		return true;
	}
/*
 * model performs a least-squares fit to the log score vs log survival function.
 * survival must be called prior to using model
 */
	bool model()
	{
		survival();
		m_fA0 = (float)3.5;
		m_fA1 = (float)-0.18;
/*
 * use default values if the statistics in the survival function is too meager
 */
		if(m_lLength == 0 || m_vlSurvive[0] < 200)	{
			return false;
		}
		float *pfX = new float[m_lLength];
		float *pfT = new float[m_lLength];
		long a = 1;
		long lMaxLimit = (long)(0.5 + m_vlSurvive[0]/2.0);
		const long lMinLimit = 10;
/*
 * find non zero points to use for the fit
 */
		a = 0;
		while(a < m_lLength && m_vlSurvive[a] > lMaxLimit)	{
			a++;
		}
		long b = 0;
		while(a < m_lLength-1 && m_vlSurvive[a] > lMinLimit)	{
			pfX[b] = (float)a;
			pfT[b] = (float)log10((double)m_vlSurvive[a]);
			b++;
			a++;
		}
/*
 * calculate the fit
 */
		double dSumX = 0.0;
		double dSumT = 0.0;
		double dSumXX = 0.0;
		double dSumXT = 0.0;
		long iMaxValue = 0;
		double dMaxT = 0.0;
		long iValues = b;
		a = 0;
		while(a < iValues)	{
			if(pfT[a] > dMaxT)	{
				dMaxT = pfT[a];
				iMaxValue = a;
			}
			a++;
		}
		a = iMaxValue;
		while(a < iValues)	{
			dSumX += pfX[a];
			dSumXX += pfX[a]*pfX[a];
			dSumXT += pfX[a]*pfT[a];
			dSumT += pfT[a];
			a++;
		}
		iValues -= iMaxValue;
		double dDelta = (double)iValues*dSumXX - dSumX*dSumX;
		if(dDelta == 0.0)	{
			delete pfX;
			delete pfT;
			return false;
		}
		m_fA0 = (float)((dSumXX*dSumT -dSumX*dSumXT)/dDelta);
		m_fA1 = (float)(((double)iValues*dSumXT - dSumX*dSumT)/dDelta);
		delete pfX;
		delete pfT;
		m_vlSurvive.clear();
		return true;
	}
/*
 * simple copy method using the = operator
 */
	mhistogram& operator=(const mhistogram &rhs)	{
		m_ulCount = rhs.m_ulCount;
		m_lLength = rhs.m_lLength;
		delete [] m_pList;
		if (rhs.m_pList == NULL)
			m_pList = NULL;
		else {
			size_t size = m_lLength * sizeof(unsigned short);
			m_pList = new unsigned short[size];
			memcpy(m_pList, rhs.m_pList, size);
		}
		m_fA0 = rhs.m_fA0;
		m_fA1 = rhs.m_fA1;
		m_dProteinFactor = rhs.m_dProteinFactor;
		m_lSum = rhs.m_lSum;
		return *this;
	}
	mhistogram& operator+=(const mhistogram &rhs)	{
		m_ulCount += rhs.m_ulCount;
		m_lLength = rhs.m_lLength;
		long a = 0;
		m_fA0 = rhs.m_fA0;
		m_fA1 = rhs.m_fA1;
		m_dProteinFactor = rhs.m_dProteinFactor;
		while(a < m_lLength)	{
			m_pList[a] += rhs.m_pList[a];
			a++;
		}
		m_lSum += rhs.m_lSum;
		return *this;
	}
};
/*
 * count_mhistogram stores information about the number of specific ion-types
 * discovered while modeling mspectrum object against a sequence collection.
 * model and survival are not used
 * starting in version 2004.04.01, this class is no longer a public mhistogram
 * this change was made to reduce memory usage, by eliminating the 
 * m_pSurvive array & reducing the size of the m_pList array
 */
class count_mhistogram
{
public:
	count_mhistogram(void) {m_lLength = 8;}
	virtual ~count_mhistogram(void) { }
	int m_lLength;

	virtual float convert(const float _f)	{
		return _f;
	}
	virtual int add(const long _c)	{
		if(_c < m_lLength && _c > -1)	{
			m_pList[_c]++;
		}
		else if(_c > -1) 	{
			m_pList[m_lLength-1]++;
		}
		else if(_c < 0) 	{
			m_pList[0]++;
		}
		return _c;
	}
/*
 * reset zeros the histogram array m_pList
 */
	virtual bool clear()	{
		int a = 0;
		while(a < m_lLength)	{
			m_pList[a] = 0;
			a++;
		}
		return true;
	}
/*
 * simple copy method using the = operator
 */
	count_mhistogram& operator=(const count_mhistogram &rhs)	{
		m_lLength = rhs.m_lLength;
		int a = 0;
		while(a < m_lLength)	{
			m_pList[a] = rhs.m_pList[a];
			a++;
		}
		return *this;
	}
	count_mhistogram& operator+=(const count_mhistogram &rhs)	{
		m_lLength = rhs.m_lLength;
		int a = 0;
		while(a < m_lLength)	{
			m_pList[a] += rhs.m_pList[a];
			a++;
		}
		return *this;
	}
/*
 * list returns a specified value from the stochastic distribution of the histogram
 */
	int list(const unsigned long _l)	{
		return m_pList[_l];
	}
/*
 * length returns the maximum length of the histogram
 */
	int length(void)	{
		return m_lLength;
	}
private:
	unsigned int m_pList[8]; // the histogram array
};
#endif //ifdef MHISTOGRAM_H

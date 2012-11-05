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

#ifndef MSPECTRUM_H
#define MSPECTRUM_H

// File version: 2004-01-07

/*
 * mspectrum objects store all information about individual tandem mass spectra. The parent ion
 * information and all m/z - intensity pairs are stored, along with msequence and mhistogram
 * objects that contain the information discovered while modeling a list of sequences.
 */

#include "mhistogram.h"
#include "msequence.h"

/*
 * mi objects contain the m/z - intensity information about a single peak in a spectrum
 */
class mi
{
public:
	mi(void) {
		m_fM = 0.0; 
		m_fI = 1.0;
	}
	virtual ~mi(void) { }

	float m_fM; // the m/z value
	float m_fI; // the intensity value
/*
 * a simple copy operation, using the = operator
 */
	mi& operator=(const mi &rhs)	{
		m_fM = rhs.m_fM;
		m_fI = rhs.m_fI;
		return *this;
	}
};

class mspectrum	{
public:
	mspectrum(void) {
		m_dRatio = 0.0;
		m_dMH = 0.0; 
		m_fI = 1.0; 
		m_fZ = 1.0; 
		m_tId = 0; 
		m_fScore = m_fScoreNext = 0.0; 
		m_fHyper = m_fHyperNext = 100.0;
		m_dExpect = 1000.0;
		m_dProteinExpect = 1000.0;
		m_bRepeat = false;
		m_bActive = true;
		m_uiType = 0;
		if(sizeof(size_t) == 8)	{
			m_tCurrentSequence = (size_t)0xAFFFFFFF;
		}
		else	{
			m_tCurrentSequence = (size_t)0xAFFF;
		}
		m_vseqBest.clear();
		m_vMI.clear();
		m_vMINeutral.clear();
		m_hHyper.clear(); // the histogram of hyper scores
		m_hConvolute.clear(); // the histogram of convolution scores
		m_chBCount.clear(); // the histogram of b-ion counts
		m_chYCount.clear(); // the histogram of y-ion counts
		m_mapScore.clear();
		m_mapCount.clear();
	}
	~mspectrum(void) { }
	
	size_t m_tId; // an identification number
	size_t	m_tCurrentSequence; // an identifier for the current sequence (used in mprocess)
	unsigned int m_uiType;
	float m_fScore; // the convolution score
	float m_fHyper; // the hyper score
	float m_fScoreNext; // next best convolution score
	float m_fHyperNext; // next best hyper score
	double m_dRatio;
	double m_dExpect; // the expectation value
	double m_dProteinExpect; // the expectation value for the associated protein
	double m_dMH; // the parent ion mass + a proton
	float m_fI; // the parent ion intensity (if available)
	float m_fZ; // the parent ion charge
	bool m_bRepeat; // a flag indicating that a better match for an individual peptide has already been found
	bool m_bActive; // a flag indicating that a spectrum is available for scoring
	vector<mi> m_vMI; // a vector containing the m/z - intensity information for fragment ions 
	vector<mi> m_vMINeutral; // a vector containing the m/z - intensity information for fragment ions 
	vector<msequence> m_vseqBest; // a vector containing the highest scoring msequence objects
	vector<double> m_vdStats;
	string m_strDescription;
	string m_strRt;
	mhistogram m_hHyper; // the histogram of hyper scores
	mhistogram m_hConvolute; // the histogram of convolution scores
	count_mhistogram m_chBCount; // the histogram of b-ion counts
	count_mhistogram m_chYCount; // the histogram of y-ion counts
	map<unsigned char,unsigned int> m_mapCount; // a map of the number of ions detected for each ion type
	map<unsigned char,float> m_mapScore; // a map of the convolution scores for each ion type
/*
 * a simple copy operator using the = operator
 */
	mspectrum& operator=(const mspectrum &rhs)	{
		m_vdStats = rhs.m_vdStats;
		m_uiType = rhs.m_uiType;
		m_dRatio = rhs.m_dRatio;
		m_hHyper = rhs.m_hHyper;
		m_hConvolute = rhs.m_hConvolute;
		m_chBCount = rhs.m_chBCount;
		m_chYCount = rhs.m_chYCount;
		m_mapCount = rhs.m_mapCount;
		m_mapScore = rhs.m_mapScore;
		m_vMI.clear();
		m_vMINeutral.clear();
		size_t a = 0;
		size_t tLength = rhs.m_vMI.size();
		while(a < tLength)	{
			m_vMI.push_back(rhs.m_vMI[a]);
			a++;
		}
		tLength = rhs.m_vMINeutral.size();
		a = 0;
		while(a < tLength)	{
			m_vMINeutral.push_back(rhs.m_vMINeutral[a]);
			a++;
		}
		m_dMH = rhs.m_dMH;
		m_fI = rhs.m_fI;
		m_fZ = rhs.m_fZ;
		m_tId = rhs.m_tId;
		m_fScore = rhs.m_fScore;
		m_fHyper = rhs.m_fHyper;
		m_fScoreNext = rhs.m_fScoreNext;
		m_fHyperNext = rhs.m_fHyperNext;
		m_dExpect = rhs.m_dExpect;
		m_dProteinExpect = rhs.m_dProteinExpect;
		m_bRepeat = rhs.m_bRepeat;
		m_vseqBest.clear();
		m_vseqBest = rhs.m_vseqBest;
		m_tCurrentSequence = rhs.m_tCurrentSequence;
		m_strDescription = rhs.m_strDescription;
		m_strRt = rhs.m_strRt;
		m_bActive = rhs.m_bActive;
		return *this;
	}
	mspectrum& operator*=(const mspectrum &rhs)	{
		m_hHyper = rhs.m_hHyper;
		m_hConvolute = rhs.m_hConvolute;
		m_chBCount = rhs.m_chBCount;
		m_chYCount = rhs.m_chYCount;
		m_mapCount = rhs.m_mapCount;
		m_mapScore = rhs.m_mapScore;
		m_bActive = rhs.m_bActive;
		m_dMH = rhs.m_dMH;
		m_fI = rhs.m_fI;
		m_fZ = rhs.m_fZ;
		m_tId = rhs.m_tId;
		m_uiType = rhs.m_uiType;
		m_fScore = rhs.m_fScore;
		m_dRatio = rhs.m_dRatio;
		m_fHyper = rhs.m_fHyper;
		m_fScoreNext = rhs.m_fScoreNext;
		m_fHyperNext = rhs.m_fHyperNext;
		m_dExpect = rhs.m_dExpect;
		m_dProteinExpect = rhs.m_dProteinExpect;
		m_bRepeat = rhs.m_bRepeat;
		m_vseqBest.clear();
		m_vseqBest = rhs.m_vseqBest;
		m_tCurrentSequence = rhs.m_tCurrentSequence;
		return *this;
	}
	mspectrum& operator+=(const mspectrum &rhs)	{
		m_vdStats = rhs.m_vdStats;
		m_hHyper += rhs.m_hHyper;
		m_hConvolute += rhs.m_hConvolute;
		m_chBCount += rhs.m_chBCount;
		m_chYCount += rhs.m_chYCount;
		m_mapCount = rhs.m_mapCount;
		m_mapScore = rhs.m_mapScore;
		m_uiType = rhs.m_uiType;
		size_t a = 0;
		size_t tLength = 0;
		m_bRepeat = false;
		if(rhs.m_fHyper == m_fHyper)	{
			m_dMH = rhs.m_dMH;
			m_fI = rhs.m_fI;
			m_fZ = rhs.m_fZ;
			m_tId = rhs.m_tId;
			m_fScore = rhs.m_fScore;
			m_fHyper = rhs.m_fHyper;
			m_dRatio = rhs.m_dRatio;
			m_fScoreNext = rhs.m_fScoreNext;
			m_fHyperNext = rhs.m_fHyperNext;
			a = 0;
			tLength = rhs.m_vseqBest.size();
			while(a < tLength)	{
				m_vseqBest.push_back(rhs.m_vseqBest[a]);
				a++;
			}
			m_tCurrentSequence = rhs.m_tCurrentSequence;
		}
		else if(rhs.m_fHyper > m_fHyper)	{
			m_dMH = rhs.m_dMH;
			m_fI = rhs.m_fI;
			m_fZ = rhs.m_fZ;
			m_tId = rhs.m_tId;
			m_fScore = rhs.m_fScore;
			m_fHyper = rhs.m_fHyper;
			m_dRatio = rhs.m_dRatio;
			m_fScoreNext = rhs.m_fScoreNext;
			m_fHyperNext = rhs.m_fHyperNext;
			m_vseqBest.clear();
			m_vseqBest = rhs.m_vseqBest;
			m_tCurrentSequence = rhs.m_tCurrentSequence;
		}
		m_dExpect = rhs.m_dExpect;
		m_dProteinExpect = rhs.m_dProteinExpect;
		return *this;
	}
	bool format()	{
		size_t tStart = m_strDescription.find('&');
		while(tStart != m_strDescription.npos)	{
			m_strDescription.replace(tStart,1,"&amp;");
			tStart = m_strDescription.find('&',tStart+1);
		}		
		tStart = m_strDescription.find('<');
		while(tStart != m_strDescription.npos)	{	
			m_strDescription.replace(tStart,1,"&lt;");
			tStart = m_strDescription.find('<',tStart+1);
		}		
		tStart = m_strDescription.find('>');
		while(tStart != m_strDescription.npos)	{
			m_strDescription.replace(tStart,1,"&gt;");
			tStart = m_strDescription.find('<',tStart+1);
		}		
		tStart = m_strDescription.find('\"');
		while(tStart != m_strDescription.npos)	{
			m_strDescription.replace(tStart,1,"&quot;");
			tStart = m_strDescription.find('\"',tStart+1);
		}		
		return true;
	}
	// contributed by Brian Pratt to correct a problem associated with not
	// clearing the m_vdStats vector when spectrum, use conditioning=no
       void clear_intensity_values() {  // clear intensities and associated stats
               m_vMI.clear();
               m_vdStats.clear();
       }
	   
	   bool reset() {
		m_fScore = 0.0; 
		m_fHyper = 100.0;
		m_dExpect = 1000.0;
		m_dProteinExpect = 1000.0;
		m_mapCount.clear();
		m_mapScore.clear();

		if(sizeof(size_t) == 8)	{
			m_tCurrentSequence = (size_t)0xAFFFFFFF;
		}
		else	{
			m_tCurrentSequence = (size_t)0xAFFF;
		}
		m_vseqBest.clear();
		return true;
	}
};
#endif

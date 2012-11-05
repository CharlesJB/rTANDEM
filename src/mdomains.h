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

#ifndef MDOMAINS_H
#define MDOMAINS_H


// File version: 2003-07-01
// File version: 2004-02-01
// File version: 2005-01-01

/*
 * maa records information about a modified residue
 * NOTE: mdomains.h has no corresponding .cpp file
 */

class maa
{
public:
	maa(void) { m_cRes = '\0'; m_lPos = 0; m_dMod = 0.0; m_cMut = '\0'; m_dPrompt = 0.0; m_strId.clear();}
	virtual ~maa(void) { }

	unsigned int m_lPos; // the sequence position of the residue (N-terminal = 0)
	float m_dMod; // mass of the modification
	char m_cRes; // single letter abbreviation for the amino acid
	char m_cMut; // single letter abbreviation for a discovered point mutation
	string m_strId; // character string representing an external accession number for a mutation/modification
	float m_dPrompt; // prompt loss from modification mass
/*
 * simple copy operation, using the = operator
 */
	maa& operator=(const maa &rhs)	{
		m_lPos = rhs.m_lPos;
		m_dMod = rhs.m_dMod;
		m_cRes = rhs.m_cRes;
		m_cMut = rhs.m_cMut;
		m_dPrompt = rhs.m_dPrompt;
		m_strId = rhs.m_strId;
		return *this;
	}
};
/*
 * mdomain objects contain information about an identified peptide sequence. each
 * identified peptide is refered to by a single mdomain object
 */
class mdomain
{
public:
	mdomain(void) {
		m_lS = 0;
		m_lE = 0; 
		m_lMissedCleaves = 0;
		m_fScore = 0.0; 
		m_fHyper = 0.0; 
		m_dMH = 0.0; 
		m_dDelta = 0.0; 
		m_vAa.clear();
		m_bUn = false;
	}
	virtual ~mdomain(void) { 
	}

	unsigned int m_lS; // the start position of the peptide in the protein sequence (N-terminus = 0)
	unsigned int m_lE; // the end position of the peptide in the protein sequence
	unsigned char m_lMissedCleaves; // missed cleavages
	float m_fScore; // the convolution score for the peptide
	float m_fHyper; // the hyper score for the peptide 
	double m_dMH; // the mass of the peptide + a proton
	// double m_dDelta replaces float m_fDelta, starting with version 2006.02.01
	// because of an issue with the accuracy of this value
	double m_dDelta; // the mass difference between the mass of the peptide and the measured mass
	bool m_bUn;

	vector<maa> m_vAa; // vector of modified amino acids
/*
 * simple copy operation, using the = operator
 */
	mdomain& operator=(const mdomain &rhs)	{
		m_lS = rhs.m_lS;
		m_lE = rhs.m_lE;
		m_lMissedCleaves = rhs.m_lMissedCleaves;
		m_dMH = rhs.m_dMH;
		m_dDelta = rhs.m_dDelta;
		m_fScore = rhs.m_fScore;
		m_fHyper = rhs.m_fHyper;
		size_t a = 0;
		m_vAa.clear();
		while(a < rhs.m_vAa.size())	{
			m_vAa.push_back(rhs.m_vAa[a]);
			a++;
		}
		m_bUn = rhs.m_bUn;
		return *this;
	}
	bool operator==(const mdomain &rhs)
	{	
		if(m_lS != rhs.m_lS)
			return false;
		if(m_lE != rhs.m_lE)
			return false;
		if(m_vAa.size() != rhs.m_vAa.size())
			return false;
		size_t a = 0;
		while(a < m_vAa.size())	{
			if(m_vAa[a].m_lPos != rhs.m_vAa[a].m_lPos)
				return false;
			if(m_vAa[a].m_dMod != rhs.m_vAa[a].m_dMod)
				return false;
			if(m_vAa[a].m_cRes != rhs.m_vAa[a].m_cRes)
				return false;
			if(m_vAa[a].m_cMut != rhs.m_vAa[a].m_cMut)
				return false;
			a++;
		}
		return true;
	}
};
#endif //ifdef MDOMAINS_H

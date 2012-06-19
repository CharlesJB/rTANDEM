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

#if !defined(SAXHANDLER_H)
#define SAXHANDLER_H

#include "expat.h"
#include <string>	// strcpm
#include <vector>
#include "mspectrum.h"

/**
* eXpat SAX parser wrapper.
*/
class SAXHandler
{
public:
	SAXHandler();
	virtual ~SAXHandler();
	void resetParser();
	/**
	* Receive notification of the start of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the start of
	* each element (such as allocating a new tree node or writing
	* output to a file).</p>
	*/
	virtual void startElement(const XML_Char *el, const XML_Char **attr);

	/**
	* Receive notification of the end of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the end of
	* each element (such as finalising a tree node or writing
	* output to a file).</p>
	*/
	virtual void endElement(const XML_Char *el);

	/**
	* Receive notification of character data inside an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method to take specific actions for each chunk of character data
	* (such as adding the data to a node or buffer, or printing it to
	* a file).</p>
	*/
	virtual void characters(const XML_Char *s, int len);

	/**
	* Open file and stream data to the SAX parser.  Must call
	* setFileName before calling this function.
	*/
	bool parse();

	inline void setFileName(const char* fileName)
	{
		m_strFileName = fileName;
	}

	// Helper functions
	inline bool isElement(const char *n1, const XML_Char *n2)
	{	return (strcmp(n1, n2) == 0); }

	inline bool isAttr(const char *n1, const XML_Char *n2)
	{	return (strcmp(n1, n2) == 0); }

	inline const char* getAttrValue(const char* name, const XML_Char **attr)
	{
		for (int i = 0; attr[i]; i += 2)
		{
			if (isAttr(name, attr[i]))
				return attr[i + 1];
		}

		return "";
	}

protected:

	XML_Parser m_parser;

	string  m_strFileName;
};

class mspectrumcondition;
class mscore;

class SAXSpectraHandler : public SAXHandler
{
public:
	SAXSpectraHandler(vector<mspectrum>& _pv, mspectrumcondition& _p, mscore& _m);
	virtual ~SAXSpectraHandler();
	enum	{
		I_Y =	0x01,
		I_B =	0x02,
		I_X =	0x04,
		I_A =	0x08,
		I_C =	0x10,
		I_Z =	0x20,
	} ion_type; // enum for referencing information about specific ion types.

protected:
	void pushSpectrum();	// Load current data into pvSpec, may have to guess charge
	void pushSpectrum(int charge);	// Load current data into pvSpec with specific charge
	void pushPeaks(bool bM = true, bool bI = true);	// Decode m_strData into a peak list

	inline void reset()
	{
		m_peaksCount = 0;
		m_precursorCharge = 0;
		m_precursorMz = 0;
		m_strActivation.clear();
		m_strRt.clear();
	}

private:
	int guessCharge();	// Guess the charge based on spectrum peaks
	void setDescription();	// Set specCurrent description based on current data
	void decode32(bool bM /*= true*/, bool bI /*= true*/);
	void decode64(bool bM /*= true*/, bool bI /*= true*/);
	unsigned long dtohl(uint32_t l, bool bNet);	// Convert from data to host long
	uint64_t dtohl(uint64_t l, bool bNet);
protected:
	string m_strData;	// For collecting character data.

	string m_strRt;
	string m_strActivation;
	bool m_bNetworkData;	// i.e. big endian
	bool m_bLowPrecision;	// i.e. 32-bit (v. 64-bit)
	bool m_bGaml;			// true if current file is GAML spectra
	string m_strDesc;		// description from GAML file <note label="Description">.....</note> or filename.scanid.charge
	double m_dSum;			// from <group> in GAML file sumI="8.41"
	double m_dMax;			// maxI="1.9242e+007"
	double m_dFactor;		// fI="192420"
	double m_dProton;
	int	m_scanNum;	// Scan number
	int	m_cidLevel;	// MS level
	int	m_peaksCount;	// Count of peaks in spectrum
	int	m_precursorCharge;	// Precursor charge
	double m_precursorMz;	// Precursor M/Z
	vector<float> m_vfM, m_vfI;	// Peak list vectors (masses and charges)

	vector<mspectrum>* m_pvSpec; //Doit pointer sur le m_vSpectra de mprocess
	mspectrumcondition* m_pSpecCond; //Doit pointer sur m_specCondition
	mscore* m_pScore; // the object that is used to score sequences and spectra

	size_t m_tId; // the id number of a spectrum
	mspectrum m_specCurrent; //Usage temporaire dans startElement et endElement
	long m_lLoaded; //if greater than max, output '.' during spectra loading
};

#endif              //SAXHANDLER_H

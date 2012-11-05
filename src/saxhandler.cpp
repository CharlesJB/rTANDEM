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

#include "stdafx.h"
#include "saxhandler.h"
#include "base64.h"
#include "mspectrumcondition.h"

#ifdef MSVC
#include <winsock2.h>   ///Pour windows
#endif


#ifdef GCC
#include <netinet/in.h>  ////Pour linux
#endif

#ifdef GCC4
#include <netinet/in.h>  ////Pour linux gcc v. 4
#endif

#ifdef OSX
#include <sys/types.h>  ////Pour mac
#include <machine/endian.h>
#endif

// Static callback handlers
static void startElementCallback(void *data, const XML_Char *el, const XML_Char **attr)
{
	((SAXHandler*) data)->startElement(el, attr);
}

static void endElementCallback(void *data, const XML_Char *el)
{
	((SAXHandler*) data)->endElement(el);
}

static void charactersCallback(void *data, const XML_Char *s, int len)
{
	((SAXHandler*) data)->characters(s, len);
}

SAXHandler::SAXHandler()
{
	m_parser = XML_ParserCreate(NULL);
	XML_SetUserData(m_parser, this);
	XML_SetElementHandler(m_parser, startElementCallback, endElementCallback);
	XML_SetCharacterDataHandler(m_parser, charactersCallback);
}


SAXHandler::~SAXHandler()
{
	XML_ParserFree(m_parser);
}

void SAXHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
}

void SAXHandler::endElement(const XML_Char *el)
{
}


void SAXHandler::characters(const XML_Char *s, int len)
{
}

bool SAXHandler::parse()
{
	FILE* pfIn = fopen(m_strFileName.data(), "r");
	if (pfIn == NULL)
	{
		cerr << "Failed to open input file '" << m_strFileName << "'.\n";
		return false;
	}
	char buffer[8192];
	int readBytes = 0;
	bool success = true;
	while (success && (readBytes = (int) fread(buffer, 1, sizeof(buffer), pfIn)) != 0)
	{
		success = (XML_Parse(m_parser, buffer, readBytes, false) != 0);
	}
	success = success && (XML_Parse(m_parser, buffer, 0, true) != 0);

	fclose(pfIn);

	if (!success)
	{
		XML_Error error = XML_GetErrorCode(m_parser);

		cerr << m_strFileName
			<< "(" << XML_GetCurrentLineNumber(m_parser) << ")"
			<< " : error " << (int) error << ": ";

		switch (error)
		{
			case XML_ERROR_SYNTAX:
			case XML_ERROR_INVALID_TOKEN:
			case XML_ERROR_UNCLOSED_TOKEN:
				cerr << "Syntax error parsing XML.";
				break;

			// TODO: Add more descriptive text for interesting errors.

			default:
				cerr << "XML Parsing error.";
				break;
		}
		cerr << "\n";
		return false;
	}
	return true;
}

/*******************************************************
 * SAXSpectraHandler class
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

SAXSpectraHandler::SAXSpectraHandler(vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
{
	m_pvSpec = &_vS;
	m_pSpecCond = &_sC;
	m_pScore = &_m;

	m_tId = 1;
	m_lLoaded = 0;
	m_bNetworkData = true;
	m_bLowPrecision = true;
	m_bGaml=false;
	m_dProton = 1.007276;

	reset();
}

SAXSpectraHandler::~SAXSpectraHandler()
{
}

void SAXSpectraHandler::pushSpectrum()
{
	int i=0;
	if(m_precursorCharge > 0)	// Known charge
	{
		pushSpectrum(m_precursorCharge);
	}
	else 	// Not sure about the charge state
	{
		// Guess charge state
		m_precursorCharge = guessCharge();

		if( m_precursorCharge == 1)
		{
			pushSpectrum(1);
		}
		else	// Multiple charge, most likely 2 or 3
		{
			pushSpectrum(2);
			m_tId += 100000000;
			pushSpectrum(3);
			m_tId -= 100000000;
		}
	}
}

void SAXSpectraHandler::pushPeaks(bool bM /*= true*/, bool bI /*= true*/)
{
	if(bM)
		m_vfM.clear();
	if(bI)
		m_vfI.clear();

	if(m_bGaml){
		int a=0;
		char *pLine = new char[m_strData.size()+1];
		char *pValue = NULL;
		strcpy(pLine,m_strData.c_str());
		pValue = pLine;
		if(bM){
			while(*pValue != '\0' && a < m_peaksCount)	{
				while(*pValue != '\0' && isspace(*pValue))
					pValue++;
				if(pValue == '\0')
					break;
				m_vfM.push_back((float)atof(pValue));
				while(*pValue != '\0' && !isspace(*pValue))
					pValue++;
				a++;
			}
		}
		else{
			while(*pValue != '\0' && a < m_peaksCount)	{
				while(*pValue != '\0' && isspace(*pValue))
					pValue++;
				if(pValue == '\0')
					break;
				m_vfI.push_back((float)atof(pValue));
				while(*pValue != '\0' && !isspace(*pValue))
					pValue++;
				a++;
			}
		}
		delete pLine;
	}
	else{
		if (m_bLowPrecision)	{
			decode32(bM,bI);
		}
		else	{
			decode64(bM,bI);
		}
	}
}

void SAXSpectraHandler::decode32(bool bM /*= true*/, bool bI /*= true*/)
{
// This code block was revised so that it packs floats correctly
// on both 64 and 32 bit machines, by making use of the uint32_t
// data type. -S. Wiley
	const char* pData = m_strData.data();
	size_t stringSize = m_strData.size();
	int setCount = 0;
	if(bM == true) {
		setCount++;
	}
	if(bI == true) {
		setCount++;
	}
	size_t size = m_peaksCount * setCount * sizeof(uint32_t);
	char* pDecoded = (char *) new char[size];
	memset(pDecoded, 0, size);
	if(m_peaksCount > 0) {
		// Base64 decoding
		// By comparing the size of the unpacked data and the expected size
		// an additional check of the data file integrity can be performed
		int length = b64_decode_mio( (char*) pDecoded , (char*) pData, stringSize );
		if(length != size) {
			cout << " decoded size " << length << " and required size " << (unsigned long)size << " dont match:\n";
			cout << " Cause: possible corrupted file.\n";
			exit(EXIT_FAILURE);
		}
	}
		// And byte order correction
	union udata {
		float fData;
		uint32_t iData;  
	} uData; 

	int n = 0;
	uint32_t* pDecodedInts = (uint32_t*)pDecoded; // cast to uint_32 for reading int sized chunks
	for(int i = 0; i < m_peaksCount; i++) {
		if(bM){
			uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
			m_vfM.push_back(uData.fData);
		}
		if(bI){
			uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
			m_vfI.push_back(uData.fData);
		}
	}
	// Free allocated memory
	delete[] pDecoded;
}

void SAXSpectraHandler::decode64(bool bM /*= true*/, bool bI /*= true*/)
{
// This code block was revised so that it packs floats correctly
// on both 64 and 32 bit machines, by making use of the uint32_t
// data type. -S. Wiley
	const char* pData = m_strData.data();
	size_t stringSize = m_strData.size();
	int setCount = 0;
	if(bM == true) {
		setCount++;
	}
	if(bI == true) {
		setCount++;
	}
	size_t size = m_peaksCount * setCount * sizeof(uint64_t);
	char* pDecoded = (char *) new char[size];
	memset(pDecoded, 0, size);
	if(m_peaksCount > 0) {
		// Base64 decoding
		// By comparing the size of the unpacked data and the expected size
		// an additional check of the data file integrity can be performed
		int length = b64_decode_mio( (char*) pDecoded , (char*) pData, stringSize );
		if(length != size) {
			cout << " decoded size " << length << " and required size " << (unsigned long)size << " dont match:\n";
			cout << " Cause: possible corrupted file.\n";
			exit(EXIT_FAILURE);
		}
	}
		// And byte order correction
	union udata {
		double fData;
		uint64_t iData;  
	} uData; 

	int n = 0;
	uint64_t* pDecodedInts = (uint64_t*)pDecoded; // cast to uint_64 for reading int sized chunks
	for(int i = 0; i < m_peaksCount; i++) {
		if(bM){
			uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
			m_vfM.push_back((float)uData.fData);
		}
		if(bI){
			uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
			m_vfI.push_back((float)uData.fData);
		}
	}
	// Free allocated memory
	delete[] pDecoded;
}

unsigned long SAXSpectraHandler::dtohl(uint32_t l, bool bNet)
{
	// mzData allows little-endian data format, so...
	// If it is not network (i.e. big-endian) data, reverse the byte
	// order to make it network format, and then use ntohl (network to host)
	// to get it into the host format.
	//if compiled on OSX the reverse is true
#ifdef OSX
	if (bNet)
	{
		l = (l << 24) | ((l << 8) & 0xFF0000) |
			(l >> 24) | ((l >> 8) & 0x00FF00);
	}
#else
	if (!bNet)
	{
		l = (l << 24) | ((l << 8) & 0xFF0000) |
			(l >> 24) | ((l >> 8) & 0x00FF00);
	}
#endif
	return l;
}

uint64_t SAXSpectraHandler::dtohl(uint64_t l, bool bNet)
{
	// mzData allows little-endian data format, so...
	// If it is not network (i.e. big-endian) data, reverse the byte
	// order to make it network format, and then use ntohl (network to host)
	// to get it into the host format.
	//if compiled on OSX the reverse is true
#ifdef OSX
	if (bNet)
	{
		l = (l << 56) | ((l << 40) & 0xFF000000000000LL) | ((l << 24) & 0x0000FF0000000000LL) | ((l << 8) & 0x000000FF00000000LL) |
			(l >> 56) | ((l >> 40) & 0x0000000000FF00LL) | ((l >> 24) & 0x0000000000FF0000LL) | ((l >> 8) & 0x00000000FF000000LL) ;
	}
#else
	if (!bNet)
	{
		l = (l << 56) | ((l << 40) & 0x00FF000000000000LL) | ((l << 24) & 0x0000FF0000000000LL) | ((l << 8) & 0x000000FF00000000LL) |
			(l >> 56) | ((l >> 40) & 0x000000000000FF00LL) | ((l >> 24) & 0x0000000000FF0000LL) | ((l >> 8) & 0x00000000FF000000LL) ;
	}
#endif
	return l;
}

void SAXSpectraHandler::pushSpectrum(int charge)
{
	int lLimit = 2000;

	/*
	* create a temporary mspectrum object
	*/
	m_specCurrent.m_tId = m_tId;
	m_specCurrent.m_strRt = m_strRt;
	m_specCurrent.m_uiType = 0;
	if(m_strActivation == "CID")	{
		m_specCurrent.m_uiType = I_Y|I_B;
	}
	else if(m_strActivation == "ETD")	{
		m_specCurrent.m_uiType = I_C|I_Z;
	}
	if(m_bGaml){
		m_specCurrent.m_dMH = m_precursorMz;
		m_specCurrent.m_vdStats.clear();
		m_specCurrent.m_vdStats.push_back(m_dSum);
		m_specCurrent.m_vdStats.push_back(m_dMax);
		m_specCurrent.m_vdStats.push_back(m_dFactor);
	}
	else{
		m_specCurrent.m_dMH = (m_precursorMz - m_dProton)*(float)charge + m_dProton;
	}
	m_specCurrent.m_fZ = (float) charge;

	mi miCurrent;
	m_specCurrent.clear_intensity_values();
	for(size_t i = 0; i < m_vfM.size(); i++)
	{
		miCurrent.m_fM = m_vfM[i];
		miCurrent.m_fI = m_vfI[i];
		// Only push if both mass and intensity are non-zero.
		if (miCurrent.m_fM && miCurrent.m_fI )
			m_specCurrent.m_vMI.push_back(miCurrent);
	}
	//if(m_bGaml && m_strDesc != "") old code, for mzML we want the spec id.
	if(m_strDesc != ""){
		m_specCurrent.m_strDescription = m_strDesc;
	}
	else{
		setDescription();
	}

	//Ajout du spectre a m_vSpectra
	if(m_pSpecCond->condition(m_specCurrent, *m_pScore) ){ // || true to force
		m_pvSpec->push_back(m_specCurrent);
		m_lLoaded++;
		if(m_lLoaded == lLimit)	{
			cout << ".";
			cout.flush();
			m_lLoaded = 0;
		}
	}
}

int SAXSpectraHandler::guessCharge()
{
	// All this small routine does is trying to guess the charge state of the precursor
	// from the ratio of the integrals of the intensities below and above m_precursorMz
	float intBelow = 0;
	float intTotal = 0;

	size_t length = m_vfM.size();
	for(size_t i = 0 ; i < length ; i++)
	{
		intTotal += m_vfI[i];
		if(m_vfM[i] < m_precursorMz)
			intBelow += m_vfI[i];   
	}

	// There is no particular reason for the 0.95. It's there just
	// to compensate for the noise.... 
	if(intTotal == 0.0 || intBelow/intTotal > 0.95)
		return 1;
	else
		return 2;
}

void SAXSpectraHandler::setDescription()
{
	m_specCurrent.m_strDescription.clear();
	size_t iPos = 0;
	size_t iSlash = m_strFileName.rfind('/');
	if (iSlash != string::npos && iSlash > iPos)
		iPos = iSlash + 1;
	size_t iBack = m_strFileName.rfind('\\');
	if (iBack != string::npos && iBack > iPos)
		iPos = iBack + 1;
	m_specCurrent.m_strDescription += m_strFileName.substr(iPos);
	m_specCurrent.m_strDescription += " scan ";
	char buffer[20];
	sprintf(buffer, "%d", m_scanNum);
	m_specCurrent.m_strDescription += buffer;
	m_specCurrent.m_strDescription += " (charge ";
	sprintf(buffer, "%d", (int)m_specCurrent.m_fZ);
	m_specCurrent.m_strDescription += buffer;
	m_specCurrent.m_strDescription += ")";
}

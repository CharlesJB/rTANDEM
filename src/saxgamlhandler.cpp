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

// File version: 2005-06-01

/*
 * saxgamlhandler.cpp contains the override methods necessary for compatibility with
 * gaml mass spectrum file formats.
 */

#include "stdafx.h"
#include "saxgamlhandler.h"

SAXGamlHandler::SAXGamlHandler( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: SAXSpectraHandler(_vS, _sC, _m)
{
	m_bMZ=false;
	m_bINT=false;
	m_bInData=false;
	m_bDesc = false;
	m_bNetworkData = false;
	m_bLowPrecision = true;
	m_bGaml=true;
	m_dSum = 100.0;
	m_dMax = 100.0;
	m_dFactor = 1.0;
}

SAXGamlHandler::~SAXGamlHandler()
{
}

void SAXGamlHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
	if(isElement("note", el) && !(strcmp("Description", getAttrValue("label", attr)))){
		m_bDesc = true;
		m_strDesc="";
	}
	else if(isElement("group", el) && !(strcmp("model", getAttrValue("type", attr)))){
		string strValue;
		strValue = getAttrValue("sumI", attr);
		if(!strValue.empty()){
			m_dSum = pow(10.0,atof(strValue.c_str()));
		}
		strValue = getAttrValue("maxI", attr);
		if(!strValue.empty()){
			m_dMax = atof(strValue.c_str());
		}
		strValue = getAttrValue("fI", attr);
		if(!strValue.empty()){
			m_dFactor = atof(strValue.c_str());
		}
		m_tId = atoi(getAttrValue("id", attr));
		m_scanNum = (int)m_tId;
		m_precursorCharge = atoi(getAttrValue("z", attr));
		m_precursorMz = atof(getAttrValue("mh", attr));
	}
	else if(isElement("GAML:trace", el) && !(strcmp("tandem mass spectrum", getAttrValue("type", attr)))){
		m_bInData = true;
	}
	else if(isElement("GAML:Xdata", el) && m_bInData){
		m_bMZ = true;
	}
	else if(isElement("GAML:Ydata", el) && m_bInData){
		m_bINT = true;
	}
	else if(isElement("GAML:values", el) && m_bInData){
		m_peaksCount = atoi(getAttrValue("numvalues", attr));
	}
}

void SAXGamlHandler::endElement(const XML_Char *el)
{
	if(isElement("note", el) && m_bDesc){
		m_bDesc=false;
	}
	else if(isElement("GAML:Xdata", el) && m_bInData){
		pushPeaks(m_bMZ, m_bINT);
		m_strData.clear();
		m_bMZ = false;
	}
	else if(isElement("GAML:Ydata", el) && m_bInData){
		pushPeaks(m_bMZ, m_bINT);
		m_strData.clear();
		m_bINT = false;
	}
	else if(isElement("GAML:trace", el) && m_bInData){
		pushSpectrum();
		reset();
		m_bInData=false;
	}
}

void SAXGamlHandler::characters(const XML_Char *s, int len)
{
	if(m_bDesc){
		m_strDesc.append(s, len);
	}
	else if((m_bMZ || m_bINT) && m_bInData){
		m_strData.append(s, len);
	}
}

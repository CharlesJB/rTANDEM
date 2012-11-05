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
// File version: 2004-03-01

/*
 * msequenceserver takes a sequence list file (in FASTA format) and uses it to load an
 * msequencecontainer object with a set of amino acid sequences and descriptions. the container
 * can be loaded repeatedly until the end-of-file is reached. the next list file path in a 
 * deque of file names is then extracted and the process continues until the last
 * sequence from the last file in the deque is read into the msequencecontainer.
 */

#include "stdafx.h"
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "xmltaxonomy.h"
#include <ctime>

p3msequenceServer::p3msequenceServer(void)
{
}

p3msequenceServer::~p3msequenceServer(void)
{
}

long p3msequenceServer::u_load_file(const string &_p,const string &_t)
{
	m_strTaxonPath = _p;
	m_strTaxon = _t;
	string strType = "peptide";
	XmlTaxonomy xmlTax;
	if(!xmlTax.load(m_strTaxonPath,m_strTaxon,strType))
		return 1;
	size_t a = 0;
	ifstream ifTest;
	m_vstrFasta.clear();
	string strTest;
	string strIn = "u_";
	while(a < xmlTax.m_vstrPaths.size())	{
		strTest = xmlTax.m_vstrPaths[a];
		strTest.insert(strTest.rfind('/')+1,strIn);
		ifTest.open(strTest.c_str());
		if(!ifTest.fail())	{
			m_dstrFasta.push_back(strTest);
			m_vstrFasta.push_back(strTest);
			ifTest.close();
		}
		ifTest.clear();
		a++;
	}
	if(m_dstrFasta.empty())	{
		return 2;
	}
	return 0;
}

unsigned long p3msequenceServer::u_maps(DESMAP &_desc,vector<msequence> &_v)
{
	unsigned long lCount = 1;
	unsigned long iLength = 0;
	msequence seqTemp;
	register char cValue = '\0';
	char *pValue = NULL;
	char *pEol = NULL;
	unsigned long lLength = 0;
	seqTemp.m_strDes = " ";
	seqTemp.m_strSeq = " ";
	m_vstrPaths.clear();
	size_t sT = 0;
	while(!m_dstrFasta.empty())	{
		m_strPath = m_dstrFasta.front();
		m_dstrFasta.pop_front();
		m_vstrPaths.push_back(m_strPath);
	/*
	* open the file	
	*/
		m_pInput = fopen(m_strPath.c_str(),"rb");
		sT = fread(m_pLine,256,1,m_pInput);
		if(strstr(m_pLine,"xbang-pro-fasta-format") != NULL)	{
			m_lFileType = XBANG;
		}
		while(!feof(m_pInput))	{
			sT = fread(&lLength,4,1,m_pInput);
#ifdef OSX
			lLength = mac_rev(lLength);
#endif
			sT = fread(m_pLine,lLength,1,m_pInput);
			seqTemp.m_strDes = m_pLine;
			if(_desc.find(seqTemp.m_strDes) != _desc.end())	{
				sT = fread(&lLength,4,1,m_pInput);
#ifdef OSX
				lLength = mac_rev(lLength);
#endif
				sT = fread(m_pLine,lLength,1,m_pInput);
				seqTemp.m_strSeq = m_pLine;
				seqTemp.m_siPath = (short int)(m_vstrPaths.size()-1);
				seqTemp.m_mapMods.clear();
				seqTemp.m_tUid = lCount;
				_v.push_back(seqTemp);
			}
			lCount++;
			iLength++;
		}
		fclose(m_pInput);
	}
	return lCount;
}

long p3msequenceServer::load_file(const string &_p,const string &_t)
{
	m_strTaxonPath = _p;
	m_strTaxon = _t;
	string strType = "peptide";
	XmlTaxonomy xmlTax;
	if(!xmlTax.load(m_strTaxonPath,m_strTaxon,strType))
		return 1;
	size_t a = 0;
	ifstream ifTest;
	m_vstrFasta.clear();
	m_dstrFasta.clear();
	while(a < xmlTax.m_vstrPaths.size())	{
		ifTest.open(xmlTax.m_vstrPaths[a].c_str());
		if(!ifTest.fail())	{
			m_dstrFasta.push_back(xmlTax.m_vstrPaths[a]);
			m_vstrFasta.push_back(xmlTax.m_vstrPaths[a]);
			ifTest.close();
		}
		ifTest.clear();
		a++;
	}
	if(m_dstrFasta.empty())	{
		return 2;
	}
	return 0;
}

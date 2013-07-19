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

// File version: 2004-02-01
// File version: 2005-01-01
/*
 * the mreport object provides the functionality to export all of the information collected
 * in mprocess to an XML file. The file path was defined in the original XML input parameter
 * file. The file produced by mreport can be used as an input parameter file to repeat
 * a protein modeling session at a later time.
 */

#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
#include "saxsaphandler.h"
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"

mreport::mreport(mscore& score)
	: m_Score(score)
{
	m_lHistogramColumns = 30;
	m_bCompress = false;
}

mreport::~mreport(void)
{
}
/*
 * Sets the value for the number of values in a GAML value list before a return
 */
bool mreport::set_columns(const long _v)
{
	if(_v < 1)
		return false;
	m_lHistogramColumns = _v;
	return true;
}
/*
 * Sets the value for the m_bCompress flag. When set, only one copy of a protein's
 * sequence is stored in the output XML file.
 */
bool mreport::set_compression(const bool _b)
{
	m_bCompress = _b;
	return m_bCompress;
}
/*
 * Checks a protein uid to see if it has already been used, storing the uid in m_vtProteins if it has not
 * been stored. This function is used to reduce the number of protein sequences stored in an XML file
 */
bool mreport::check_proteins(const size_t _t)
{
	if(!m_bCompress)
		return true;
	if(m_setProteins.find(_t) == m_setProteins.end())	{
		m_setProteins.insert(_t);
		return true;
	}
	return false;
}
/*
 * end finishes the XML document
 */
bool mreport::end(void)
{
	if(m_ofOut.fail())	{
		return false;
	}
	m_ofOut << "</bioml>\n";
	m_ofOut.close();
	return true;
}
/*
 * endgroup finishes the group object associated with a particular mspectrum object
 */
bool mreport::endgroup()
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	m_ofOut << "</group>\n";
	return true;
}
/*
 * group starts the group associated with outputing the information from a particular
 * mspectrum object
 */
bool mreport::group(const mspectrum &_s)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	char *pLine = new char[256];
	string strTemp;
	size_t tId = _s.m_tId;
	while(tId > 100000000)	{
		tId -= 100000000;
	}
	if(_s.m_vseqBest.size() == 0)	{
		m_ofOut << "<group id=\"" << (unsigned long)tId << "\" ";
		sprintf(pLine,"%.6lf",_s.m_dMH);
		m_ofOut << "mh=\"" << pLine << "\" ";
		m_ofOut << "z=\"" << (long)_s.m_fZ << "\" ";
		m_ofOut << "rt=\"" << _s.m_strRt.c_str() << "\" ";
		m_ofOut << "expect=\"1000\" ";
		sprintf(pLine,"%.2lf",log10(_s.m_vdStats[0]));
		m_ofOut << "label=\"" << "no model obtained" << "\" type=\"model\" ";
		m_ofOut << "sumI=\"" << pLine << "\" maxI=\"" << _s.m_vdStats[1] << "\" fI=\"" << _s.m_vdStats[2] << "\" ";
		m_ofOut << "act=\"" << _s.m_uiType << "\" ";
		m_ofOut << " >\n";
	}
	else	{
		m_ofOut << "<group id=\"" << (unsigned long)tId << "\" ";
		sprintf(pLine,"%.6lf",_s.m_dMH);
		m_ofOut << "mh=\"" << pLine << "\" ";
		m_ofOut << "z=\"" << (long)_s.m_fZ << "\" ";
		m_ofOut << "rt=\"" << _s.m_strRt.c_str() << "\" ";
		sprintf(pLine,"%.1e",_s.m_dExpect);
		m_ofOut << "expect=\"" << pLine << "\" ";
		string strV = _s.m_vseqBest[0].m_strDes;
		format_text(strV);
		get_short_label(strV,pLine,80,255);
		m_ofOut << "label=\"" << pLine << "\" type=\"model\" ";
		sprintf(pLine,"%.2lf",log10(_s.m_vdStats[0]));
		m_ofOut << "sumI=\"" << pLine << "\" maxI=\"" << _s.m_vdStats[1] << "\" fI=\"" << _s.m_vdStats[2] << "\" ";
		m_ofOut << "act=\"" << _s.m_uiType << "\" ";
		m_ofOut << ">\n";
	}
	delete pLine;
	return true;
}

bool mreport::get_short_label(const string &_s,char *_p,const unsigned long _l,const unsigned long _t)
{
	size_t lLength = _s.size();
	size_t a = 0;
	size_t tMax = _t - 5;
	while(a < lLength && a < _l && strchr("\r\n",_s[a]) == NULL)	{
		_p[a] = _s[a];
		a++;
	}
	if(strchr("\r\n",_s[a]))	{
		_p[a] = '\0';
		return true;
	}
	while(a < lLength && !isspace(_s[a]) && a < tMax)	{
		_p[a] = _s[a];
		a++;
	}
	if(a != lLength)	{
//		cout << _s.c_str() << "||";
		// Rprintf("%s||", _s.c_str());
		//cout.flush();
		_p[a] = '.';
		a++;
		_p[a] = '.';
		a++;
		_p[a] = '.';
		a++;
	}
	_p[a] = '\0';
	return true;
}
/*
 * histogram outputs all of the histrograms contained in an mspectrum object. GAML
 * notation is used to output the histograms
 */
bool mreport::histogram(mspectrum &_s)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	__int64_t tId = _s.m_tId;
	while(tId > 100000000)	{
		tId -= 100000000;
	}
	m_ofOut << "<group label=\"supporting data\" type=\"support\">\n";
	m_ofOut << "<GAML:trace label=\"" << tId << ".hyper\" type=\"hyperscore expectation function\">\n";
	m_ofOut << "<GAML:attribute type=\"a0\">" << _s.m_hHyper.a0() << "</GAML:attribute>\n";
	m_ofOut << "<GAML:attribute type=\"a1\">" << _s.m_hHyper.a1() << "</GAML:attribute>\n";
	m_ofOut << "<GAML:Xdata label=\"" << tId << ".hyper\" units=\"score\">\n";
	long a = _s.m_hHyper.length()-1;
	_s.m_hHyper.survival();
	while(a >= 0 && _s.m_hHyper.survive(a) < 1)	{
		a--;
	}
	if(a == _s.m_hHyper.length()-1)
		a = _s.m_hHyper.length()-2;
	long lLength = a+2;
	long lCount = 0;
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	while(a < lLength)	{
		m_ofOut << a;
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Xdata>\n";
	m_ofOut << "<GAML:Ydata label=\"" << tId << ".hyper\" units=\"counts\">\n";
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << _s.m_hHyper.survive(a);
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	_s.m_hHyper.clear_survive();
	m_ofOut << "\n</GAML:values>\n</GAML:Ydata>\n</GAML:trace>\n";
	m_ofOut << "<GAML:trace label=\"" << tId << ".convolute\" type=\"convolution survival function\">\n";
	m_ofOut << "<GAML:Xdata label=\"" << tId << ".convolute\" units=\"score\">\n";
	a = _s.m_hConvolute.length()-1;
	_s.m_hConvolute.survival();
	while(a >= 0 && _s.m_hConvolute.survive(a) < 1)	{
		a--;
	}
	if(a == _s.m_hConvolute.length()-1)
		a = _s.m_hConvolute.length()-2;
	lLength = a+2;
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << a;
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Xdata>\n";
	m_ofOut << "<GAML:Ydata label=\"" << tId << ".convolute\" units=\"counts\">\n";
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << _s.m_hConvolute.survive(a);
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	_s.m_hConvolute.clear_survive();
	m_ofOut << "\n</GAML:values>\n</GAML:Ydata>\n</GAML:trace>\n";
	m_ofOut << "<GAML:trace label=\"" << tId << ".b\" type=\"b ion histogram\">\n";
	m_ofOut << "<GAML:Xdata label=\"" << tId << ".b\" units=\"number of ions\">\n";
	a = _s.m_chBCount.length()-1;
	while(a >= 0 && _s.m_chBCount.list(a) < 1)	{
		a--;
	}
	if(a == _s.m_chBCount.length()-1)
		a = _s.m_chBCount.length()-2;
	lLength = a+2;
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << a;
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Xdata>\n";
	m_ofOut << "<GAML:Ydata label=\"" << tId << ".b\" units=\"counts\">\n";
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << _s.m_chBCount.list(a);
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Ydata>\n</GAML:trace>\n";
	m_ofOut << "<GAML:trace label=\"" << tId << ".y\" type=\"y ion histogram\">\n";
	m_ofOut << "<GAML:Xdata label=\"" << tId << ".y\" units=\"number of ions\">\n";
	a = _s.m_chYCount.length()-1;
	while(a >= 0 && _s.m_chYCount.list(a) < 1)	{
		a--;
	}
	if(a == _s.m_chYCount.length()-1)
		a = _s.m_chYCount.length()-2;
	lLength = a+2;
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << a;
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Xdata>\n";
	m_ofOut << "<GAML:Ydata label=\"" << tId << ".y\" units=\"counts\">\n";
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< lLength << "\">\n";
	a = 0;
	lCount = 0;
	while(a < lLength)	{
		m_ofOut << _s.m_chYCount.list(a);
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Ydata>\n</GAML:trace>\n";
	m_ofOut << "\n</group>\n";
	return true;
}
/*
 * info outputs the input parameters derived from the input file and the default input file
 * (if used). This record should be sufficient to repeat this protein modeling session exactly, without reference to
 * the original input files.
 * Two group nodes are used to record this information. The "input parameters" group contains all
 * of the input nodes that were used to create this output file. The "unused input parameters" group
 * contains all of the input parameters that were either unsupported, mis-entered or otherwise
 * ignored.
 */

bool mreport::info(XmlParameter &_x)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	pair<string,string> pairValue;
	xMap::iterator itStart = _x.m_mapParam.begin();
	xMap::iterator itEnd = _x.m_mapParam.end();
	size_t a = 0;
	string strValue;
	m_ofOut << "<group label=\"input parameters\" type=\"parameters\">\n";
	size_t tUnused = 0;
	while(itStart != itEnd)	{
		pairValue = *itStart;
		if(_x.m_mapUsed[pairValue.first])	{
			strValue = pairValue.second;
			m_ofOut << "\t<note type=\"input\" label=\"" << pairValue.first.c_str() << "\"";
			m_ofOut << ">";
			a = 0;
			while(a < strValue.size())	{
				if(strValue[a] == '<')	{
					m_ofOut << "&lt;";
				}
				else if(strValue[a] == '>')	{
					m_ofOut << "&gt;";
				}
				else if(strValue[a] == '\"')	{
					m_ofOut << "&quot;";
				}
				else	{
					m_ofOut << strValue[a];
				}
				a++;
			}
			m_ofOut << "</note>\n";
		}
		else	{
			tUnused++;
		}
		itStart++;
	}
	m_ofOut << "</group>\n";
	if(tUnused == 0)
		return true;
	itStart = _x.m_mapParam.begin();
	itEnd = _x.m_mapParam.end();
	a = 0;
	m_ofOut << "<group label=\"unused input parameters\"  type=\"parameters\">\n";
	while(itStart != itEnd)	{
		pairValue = *itStart;
		if(!_x.m_mapUsed[pairValue.first])	{
			strValue = pairValue.second;
			m_ofOut << "\t<note type=\"input\" label=\"" << pairValue.first.c_str() << "\"";
			m_ofOut << ">";
			a = 0;
			while(a < strValue.size())	{
				if(strValue[a] == '<')	{
					m_ofOut << "&lt;";
				}
				else if(strValue[a] == '>')	{
					m_ofOut << "&gt;";
				}
				else if(strValue[a] == '\"')	{
					m_ofOut << "&quot;";
				}
				else	{
					m_ofOut << strValue[a];
				}
				a++;
			}
			m_ofOut << "</note>\n";
		}
		itStart++;
	}
	m_ofOut << "</group>\n";
	return true;
}
/*
 * performance is called to output a group node that contains any
 * performance parameters from the protein modeling session that might be of interest
 * after the protein modeling session. Any number of additional output parameters can
 * be included here. An XmlParameter object is used to supply the information
 * required to create a series of note objects that contain the informtion.
 */
bool mreport::performance(XmlParameter &_x)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	pair<string,string> pairValue;
	xMap::iterator itStart = _x.m_mapParam.begin();
	xMap::iterator itEnd = _x.m_mapParam.end();
	size_t a = 0;
/*
 * create the group object
 */
	string strValue;
	m_ofOut << "<group label=\"performance parameters\" type=\"parameters\">\n";
	while(itStart != itEnd)	{
		pairValue = *itStart;
		strValue = pairValue.second;
/*
 * create the individual notes
 */
		m_ofOut << "\t<note label=\"" << pairValue.first << "\">";
		a = 0;
		while(a < strValue.size())	{
			if(strValue[a] == '<')	{
				m_ofOut << "&lt;";
			}
			else if(strValue[a] == '>')	{
				m_ofOut << "&gt;";
			}
			else if(strValue[a] == '\"')	{
				m_ofOut << "&quot;";
			}
			else	{
				m_ofOut << strValue[a];
			}
			a++;
		}
		m_ofOut << "</note>\n";
		itStart++;
	}
/*
 * finish the group
 */
	m_ofOut << "</group>\n";
	return true;
}
/*
 * write masses to the results file if residue masses are not default values
*/
bool mreport::masses(msequtilities &_p)
{
	if(!_p.is_modified())	{
		return false;
	}
	char cAa = 'A';
/*
 * create the group object
 */
	string strValue;
	char *pLine = new char[256];
	m_ofOut << "<group label=\"residue mass parameters\" type=\"parameters\">\n";
	while(cAa <= 'Z')	{
		int index = cAa;
		sprintf(pLine,"\t<aa type=\"%c\" mass=\"%.6lf\" />\n",cAa,_p.m_pdAaMass[index]);
		m_ofOut << pLine;
		cAa++;
	}
	sprintf(pLine,"\t<molecule type=\"NH3\" mass=\"%.6lf\" />\n",_p.m_dAmmonia);
	m_ofOut << pLine;
	sprintf(pLine,"\t<molecule type=\"H2O\" mass=\"%.6lf\" />\n",_p.m_dWater);
	m_ofOut << pLine;
/*
 * finish the group
 */
	m_ofOut << "</group>\n";
	delete pLine;
	return true;
}
/*
 * spectrum is used to output a complete spectrum, using GAML notation.
 */
bool mreport::spectrum(mspectrum &_s,string &_f)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	__int64_t tId = _s.m_tId;
	while(tId > 100000000)	{
		tId -= 100000000;
	}
	m_ofOut << "<group type=\"support\" label=\"fragment ion mass spectrum\">\n";
	if(!_f.empty())	{
		m_ofOut << "<file type=\"spectra\" URL=\"" << _f.c_str() << "\" />\n";
	}
	if(!_s.m_strDescription.empty())	{
		_s.format();
		m_ofOut << "<note label=\"Description\">" << _s.m_strDescription.c_str() << "</note>\n";
	}
	m_ofOut << "<GAML:trace id=\"" << tId << "\" label=\"" << tId << ".spectrum\" type=\"tandem mass spectrum\">\n";
	m_ofOut << "<GAML:attribute type=\"M+H\">" << _s.m_dMH << "</GAML:attribute>\n";
	m_ofOut << "<GAML:attribute type=\"charge\">" << _s.m_fZ << "</GAML:attribute>\n";
	m_ofOut << "<GAML:Xdata label=\"" << tId << ".spectrum\" units=\"MASSTOCHARGERATIO\">\n";
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< (unsigned long)_s.m_vMI.size() << "\">\n";
	size_t a = 0;
	const size_t tLength = _s.m_vMI.size();
	long lCount = 0;
	while(a < tLength)	{
		m_ofOut << _s.m_vMI[a].m_fM;
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Xdata>\n";
	m_ofOut << "<GAML:Ydata label=\"" << tId << ".spectrum\" units=\"UNKNOWN\">\n";
	m_ofOut << "<GAML:values byteorder=\"INTEL\" format=\"ASCII\" numvalues=\"" 
		<< (unsigned long)_s.m_vMI.size() << "\">\n";
	a = 0;
	lCount = 0;
	char *pLine = new char[256];
	while(a < tLength)	{
		sprintf(pLine,"%.0f",_s.m_vMI[a].m_fI);
		m_ofOut << pLine;
		lCount++;
		if(lCount == m_lHistogramColumns)	{
			m_ofOut << "\n";
			lCount = 0;
		}
		else	{
			m_ofOut  << " ";
		}
		a++;
	}
	m_ofOut << "\n</GAML:values>\n</GAML:Ydata>\n</GAML:trace>\n</group>";
	delete pLine;
	return true;
}
/*
 * sequence outputs the information about an identified protein sequence, using
 * bioml notation
 */
bool mreport::sequence(mspectrum &_s,const bool _b,vector<string> &_p,map<string,string> &_ann)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	size_t tRow = 0;
	size_t tSpace = 0;
	char *pLine = new char[256];
	__int64_t tId = _s.m_tId;
	while(tId > 100000000)	{
		tId -= 100000000;
	}
	string strTemp;
	map<string,string>::iterator itMod;
	while(a < _s.m_vseqBest.size())	{
		m_ofOut << "<protein expect=\"";
		sprintf(pLine,"%.1lf",_s.m_vseqBest[a].m_dExpect);
		if(strstr(pLine,"-1.$"))	{
				strcpy(pLine,"-5999.0");
		}
		m_ofOut << pLine << "\"";
		m_ofOut << " id=\"" << tId << "." << (unsigned long)(a+1) << "\"";
		m_ofOut << " uid=\"" << (unsigned long)_s.m_vseqBest[a].m_tUid << "\" ";
		string strV = _s.m_vseqBest[a].m_strDes;
		format_text(strV);
		get_short_label(strV,pLine,80,250);

		m_ofOut << "label=\"" << pLine << "\" ";
		sprintf(pLine,"sumI=\"%.2lf\" ",log10(_s.m_vseqBest[a].m_fIntensity));
		m_ofOut << pLine;
		itMod = _ann.find(_s.m_vseqBest[a].m_strDes);
		if(itMod != _ann.end())	{
			sprintf(pLine," annotation=\"%s\"",(itMod->second).c_str());
			m_ofOut << pLine;
		}
		m_ofOut <<  ">\n";
		m_ofOut << "<note label=\"description\">";
		m_ofOut << strV.c_str() << "</note>\n";
		m_ofOut << "<file type=\"peptide\" URL=\"" << _p[_s.m_vseqBest[a].m_siPath] << "\"/>\n";
		m_ofOut << "<peptide start=\"1\" end=\"" << (unsigned long)_s.m_vseqBest[a].m_strSeq.size() << "\">\n";
		c = 0;
		tRow = 0;
		tSpace = 0;
		m_ofOut << "\t";
		if(m_maptUids.find(_s.m_vseqBest[a].m_tUid) == m_maptUids.end())	{
			m_maptUids[_s.m_vseqBest[a].m_tUid] = true;
		}
		if(check_proteins(_s.m_vseqBest[a].m_tUid))	{
			while(_b && c < _s.m_vseqBest[a].m_strSeq.size())	{
				m_ofOut << _s.m_vseqBest[a].m_strSeq[c];
				tRow++;
				tSpace++;
				if(tRow == 50)	{
					tSpace = 0;
					tRow = 0;
					m_ofOut << "\n\t";
				}
				else if(tSpace == 10)	{
					tSpace = 0;
					m_ofOut << " ";
				}
				c++;
			}
			m_ofOut << "\n";
		}
		b = 0;
		double dDelta;
		while(b <  _s.m_vseqBest[a].m_vDomains.size())	{
			m_ofOut << "<domain id=\"" << tId << "." << (unsigned long)(a+1) << "." << (unsigned long)(b+1) << "\" ";
			m_ofOut << "start=\"" << (unsigned long)(_s.m_vseqBest[a].m_vDomains[b].m_lS+1) << "\" ";
			m_ofOut << "end=\"" << (unsigned long)(_s.m_vseqBest[a].m_vDomains[b].m_lE+1) << "\" ";
			sprintf(pLine,"%.1e",_s.m_hHyper.expect(m_Score.hconvert(_s.m_vseqBest[a].m_vDomains[b].m_fHyper)));
			m_ofOut << "expect=\"" << pLine << "\" ";
			// dDelta was subsituted for _s.m_vseqBest[a].m_fDelta in version 2006.02.01 to
			// improve the accuracy of the mass difference calculation
			dDelta = _s.m_dMH - _s.m_vseqBest[a].m_vDomains[b].m_dMH;
			if(fabs(dDelta) > 0.01)	{
				sprintf(pLine,"%.3lf",_s.m_vseqBest[a].m_vDomains[b].m_dMH);
			}
			else	{
				sprintf(pLine,"%.4lf",_s.m_vseqBest[a].m_vDomains[b].m_dMH);
			}
			m_ofOut << "mh=\"" << pLine << "\" ";
			if(fabs(dDelta) > 0.01)	{
				sprintf(pLine,"%.3lf",dDelta);
			}
			else	{
				sprintf(pLine,"%.4lf",dDelta);
			}
			m_ofOut << "delta=\"" << pLine << "\" ";
			m_Score.report_score(pLine, _s.m_vseqBest[a].m_vDomains[b].m_fHyper);
			m_ofOut << "hyperscore=\"" << pLine << "\" ";
			m_Score.report_score(pLine, _s.m_fHyperNext);
			m_ofOut << "nextscore=\"" << pLine << "\" ";
			m_Score.report_score(pLine, _s.m_mapScore[mscore::T_Y]+_s.m_mapScore[mscore::T_Z]);
			m_ofOut << "y_score=\"" << pLine << "\" ";
			m_ofOut << "y_ions=\"" << _s.m_mapCount[mscore::T_Y]+_s.m_mapCount[mscore::T_Z] << "\" ";
			m_Score.report_score(pLine, _s.m_mapScore[mscore::T_B]+_s.m_mapScore[mscore::T_C]);
			m_ofOut << "b_score=\"" << pLine << "\" ";
			m_ofOut << "b_ions=\"" << _s.m_mapCount[mscore::T_B] + _s.m_mapCount[mscore::T_C] << "\" ";
			get_pre(_s.m_vseqBest[a].m_strSeq,strTemp,_s.m_vseqBest[a].m_vDomains[b].m_lS,_s.m_vseqBest[a].m_vDomains[b].m_lE);
			m_ofOut << "pre=\"" << strTemp.c_str()  << "\" ";
			get_post(_s.m_vseqBest[a].m_strSeq,strTemp,_s.m_vseqBest[a].m_vDomains[b].m_lS,_s.m_vseqBest[a].m_vDomains[b].m_lE);
			m_ofOut << "post=\"" << strTemp.c_str() << "\" ";
			m_ofOut << "seq=\"" << _s.m_vseqBest[a].m_strSeq.substr(_s.m_vseqBest[a].m_vDomains[b].m_lS,
									_s.m_vseqBest[a].m_vDomains[b].m_lE-_s.m_vseqBest[a].m_vDomains[b].m_lS+1).c_str() << "\" ";
			m_ofOut << "missed_cleavages=\"" << (int)_s.m_vseqBest[a].m_vDomains[b].m_lMissedCleaves
									<< "\">\n";
			c = 0;
			while(c < _s.m_vseqBest[a].m_vDomains[b].m_vAa.size())	{
				m_ofOut << "<aa type=\"" << _s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_cRes << "\" at=\"" << _s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_lPos+1 << "\" ";
				sprintf(pLine,"%.5lf",_s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_dMod);
				m_ofOut << "modified=\"" << pLine << "\" ";
				if(_s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_cMut != '\0' && _s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_cMut != _s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_cRes)	{
					m_ofOut << "pm=\"" << _s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_cMut << "\" ";
				}
				if(_s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_strId.size() != 0)	{
					m_ofOut << "id=\"" << _s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_strId.c_str() << "\" ";
				}
				if(fabs(_s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_dPrompt) > 0.1)	{
					sprintf(pLine,"%.5lf",_s.m_vseqBest[a].m_vDomains[b].m_vAa[c].m_dPrompt);
					m_ofOut << "prompt=\"" << pLine << "\" ";
				}
				m_ofOut << "/>\n";
				c++;
			}
			m_ofOut << "</domain>\n";
			b++;
		}
		m_ofOut << "</peptide>\n</protein>\n";
		a++;
	}
	delete pLine;
	return true;
}

bool mreport::get_pre(const string &_s,string &_p,const size_t _b,const size_t _e)
{
	long lStart = (long)_b;
	lStart -= 4;
	_p.erase(_p.begin(),_p.end());
	if(lStart < 0)	{
		lStart = 0;
		_p = '[';
	}
	while(lStart < (long)_b)	{
		_p += _s[lStart];
		lStart++;
	}
	return true;
}

bool mreport::get_post(const string &_s,string &_p,const size_t _b,const size_t _e)
{
	size_t tEnd = (long)_e;
	tEnd += 5;
	_p.erase(_p.begin(),_p.end());
	if(tEnd >= _s.size())	{
		tEnd = _s.size();
	}
	size_t a = _e + 1;
	while(a < tEnd)	{
		_p += _s[a];
		a++;;
	}
	if(a == _s.size())	{
		_p += ']';
	}
	return true;
}

/*
 * start is called first to set up the XML header
 */
bool mreport::start(XmlParameter &_x)
{
	string strKey = "output, path";
	string strValue;
	_x.get(strKey,strValue);
	if(strValue.size() == 0)
		return false;
	string strPath = strValue;
	strKey = "output, path hashing";
	_x.get(strKey,strValue);
	if(strValue == "yes")	{
		size_t tStart = strPath.rfind('.');
		if(tStart != strPath.npos)	{
			tStart++;
		    time_t tValue;
			time(&tValue);
			char pLine[256];
			struct tm *ptmValue = localtime(&tValue);
			if(ptmValue == NULL)	{
				ptmValue = gmtime(&tValue);
			}
			if(ptmValue != NULL)	{
				strftime(pLine, 255,"%Y_%m_%d_%H_%M_%S.t.",ptmValue);
			}
			else	{
				srand((unsigned)time(NULL));
				sprintf(pLine,"%i.t",rand());
			}
			strPath.insert(tStart,pLine);
		}
	}
	m_pathName = strPath; // rTANDEM
	m_ofOut.open(strPath.c_str());
	if(m_ofOut.fail())	{
		return false;
	}
	m_ofOut << "<?xml version=\"1.0\"?>\n";
	strKey = "output, xsl path";
	_x.get(strKey,strValue);
	if(strValue.size() > 0)	{
		m_ofOut << "<?xml-stylesheet type=\"text/xsl\" href=\"";
		m_ofOut << strValue.c_str() << "\"?>\n";
	}
	strKey = "output, title";
	_x.get(strKey,strValue);
	if(strValue.size() > 0)	{
		m_ofOut << "<bioml xmlns:GAML=\"http://www.bioml.com/gaml/\" ";
		m_ofOut << "label=\"" << strValue.c_str() << "\">\n";
	}
	else	{
		strKey = "spectrum, path";
		_x.get(strKey,strValue);
		m_ofOut << "<bioml xmlns:GAML=\"http://www.bioml.com/gaml/\" ";
		m_ofOut << "label=\"models from '" << strValue.c_str() << "'\">\n";
	}
	return true;
}

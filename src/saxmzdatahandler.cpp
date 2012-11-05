/************************************************************
 *              SAXMzdataHandler.cpp
 *
 * Premiere version janvier 2005
 * Patrick Lacasse
 * placasse@mat.ulaval.ca
 *
 * 3/11/2005 (Brendan MacLean): Use eXpat SAX parser, and create SAXSpectraHandler
 *
 * November 2005
 * Fredrik Levander 
 * Fredrik.Levander@elmat.lth.se 
 * A few changes to handle MzData 1.05.
 *
 * Updated to handle version 1.04 and 1.05. (Rob Craig)
 *
 *
 * See http://psidev.sourceforge.net/ms/#mzdata for
 * mzData schema information.
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

#include "stdafx.h"
#include "saxmzdatahandler.h"

SAXMzdataHandler::SAXMzdataHandler( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: SAXSpectraHandler(_vS, _sC, _m)
{
	m_bInmzArrayBinary=false;
	m_bInintenArrayBinary=false;
	m_bInData=false;
	m_bInMsLevel2=false;
}

SAXMzdataHandler::~SAXMzdataHandler()
{
}

void SAXMzdataHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
	// Let's see what	element	just started
	// We	use:
	//   spectrum	for	the	id number
	//   acqInstrument for the msLevel
	//   cvParam for the polarity	and	the	mz of the parent ion
	//   data	for	the	values of the spectrum

	if (isElement("spectrum", el))
	{
		m_scanNum =	atoi(getAttrValue("id", attr));
		m_tId = m_scanNum;
		while(m_sId.find(m_tId) != m_sId.end())	{
			m_tId++;
		}
		m_sId.insert(m_tId);

	}
	//mzdata v. 1.04 uses acqInstrument and mzdata v. 1.05 uses spectrumInstrument
	else if (isElement("spectrumInstrument", el) || isElement("acqInstrument", el))
	{
		//first	attribute is msLevel
		if(	(m_cidLevel	= atoi(getAttrValue("msLevel", attr))) == 2 )
		{
			m_bInMsLevel2 = true;

			reset();	// Clean up	for	the	next scan

		}
	}
	else if (isElement("cvParam", el))
	{
		const char* name = getAttrValue("name", attr);
		const char* value = getAttrValue("value", attr);
		if (!strcmp(name, "polarity")){
			if(!strcmp(value, "+"))  
				m_precursorCharge	= 0;
			else{
				m_precursorCharge = atoi(value);
				if(m_precursorCharge ==	0)
					m_precursorCharge = 2;
			}
		}
		else if (!strcmp(name, "mz") )	{
			m_precursorMz =	atof(value);
		}
		else if (!strcmp(name, "m/z") )	{
			m_precursorMz =	atof(value);
		}
		else if (!strcmp(name, "ChargeState")){
			m_precursorCharge = atoi(value);
		}
		else if (!strcmp(name, "selected ion m/z")){
			m_precursorMz =	atof(value);
		}
		else if (!strcmp(name, "Charge State")){
			m_precursorCharge = atoi(value);
		}
		else if (!strcmp(name, "Mass To Charge Ratio") ){
			m_precursorMz =	atof(value);
		}
		else if (!strcmp(name, "MassToChargeRatio") ){
			m_precursorMz =	atof(value);
		}
	}
	else if (isElement("mzArrayBinary", el)){
		m_bInmzArrayBinary = true;
		startPeakListBinary(attr);
	}
	else if (isElement("intenArrayBinary", el))	{
		m_bInintenArrayBinary = true;
		startPeakListBinary(attr);
	}
	else if (isElement("data", el))	{
		m_bInData	= true;
		m_peaksCount = atoi(getAttrValue("length", attr));
		if(strlen(getAttrValue("endian", attr)) > 0)	{
			m_bNetworkData = (strcmp("big", getAttrValue("endian", attr)) != 0);
		}
		if(strlen(getAttrValue("precision", attr)) > 0)	{
			m_bLowPrecision = (strcmp("64", getAttrValue("precision", attr)) != 0);
		}
	}
}

void SAXMzdataHandler::startPeakListBinary(const XML_Char **attr)
{
	if(strlen(getAttrValue("endian", attr)) > 0)	{
		m_bNetworkData = (strcmp("little", getAttrValue("endian", attr)) != 0);
	}
	if(strlen(getAttrValue("precision", attr)) > 0)	{
		m_bLowPrecision = (strcmp("64", getAttrValue("precision", attr)) != 0);
	}
}

void SAXMzdataHandler::endElement(const XML_Char *el)
{
	if(isElement("mzArrayBinary", el))
		m_bInmzArrayBinary = false;
	else if(isElement("intenArrayBinary", el))
		m_bInintenArrayBinary = false;
	else if(isElement("data", el))
	{
		processData();
		m_bInData = false;
	}
	else if(isElement("spectrum", el) && m_bInMsLevel2)
	{
		pushSpectrum();
		m_bInMsLevel2 = false;
	}
}

void SAXMzdataHandler::characters(const XML_Char *s, int len)
{
	if((m_bInmzArrayBinary || m_bInintenArrayBinary) && m_bInMsLevel2 && m_bInData)
	{
		m_strData.append(s, len);
	}
}

void SAXMzdataHandler::processData()
{
	if((m_bInmzArrayBinary || m_bInintenArrayBinary) && m_bInMsLevel2 && m_bInData)
	{
		pushPeaks(m_bInmzArrayBinary, m_bInintenArrayBinary);
	}

	m_strData.clear();
}

#if !defined(SAXMZMLHANDLER_H)
#define SAXMZMLHANDLER_H

#include "saxhandler.h"
#include <string>
#include <set>

/**
* Uses eXpat SAX parser to parse mzData XML.
*/
class SAXMzmlHandler : public SAXSpectraHandler
{
public:
	SAXMzmlHandler( vector<mspectrum>& _pv, mspectrumcondition& _p, mscore& _m);
	~SAXMzmlHandler();

	// -----------------------------------------------------------------------
	//  Overrides of SAXHandler functions
	// -----------------------------------------------------------------------
	/**
	* Receive notification of the start of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the start of
	* each element (such as allocating a new tree node or writing
	* output to a file).</p>
	*/
	void startElement(const XML_Char *el, const XML_Char **attr);

	/**
	* Receive notification of the end of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the end of
	* each element (such as finalising a tree node or writing
	* output to a file).</p>
	*/
	void endElement(const XML_Char *el);

	/**
	* Receive notification of character data inside an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method to take specific actions for each chunk of character data
	* (such as adding the data to a node or buffer, or printing it to
	* a file).</p>
	*/
	void characters(const XML_Char *s, int len);

private:

	/**
	* Handle attributes for a peaksListBinaryType element.
	*/
	void startPeakListBinary(const XML_Char **attr);
	set<size_t> m_sId;

	/**
	* Handle data collected in the characters() function.
	*/
	void processData();
	
	/**
	* Process cv params
	*/
	void processCVParam(const char* name, const char* accession, const char* value);

private:
	// Flags indicating parser is inside a particular tag.
	bool m_bInRefGroup;
	bool m_bInMsLevel2;
	bool m_bInmzArrayBinary;
	bool m_bInintenArrayBinary;
	bool m_bInData;
	string m_ccurrentRefGroupName;
	class cvParam
	{
	public:
		string refGroupName;
		string name;
		string accession;
		string value;
	};
	vector <cvParam> m_refGroupCvParams;

};
#endif              //SAXMZMLHANDLER_H

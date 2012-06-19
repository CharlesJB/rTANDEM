/**
 * dataLoader.h
 *
 * Author: Charles Joly Beauparlant
 * Creation Date: 18 June 2012
 *
 */

#idfndef DATALOADER_H
#define DATALOADER_H

#include <Rcpp.h>
#include <map>
#include <string>

class dataLoader {
public:
	static void convertSEXP(SEXP RData, std::map<std::string, std::string>* ptrMap);
	
private:
	dataLoader() { }
};

#endif /* DATALOADER_H */

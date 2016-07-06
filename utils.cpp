//
//  utils.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-05.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>
#include "utils.h"


void stopif(bool condition, string error_msg,
			int error_code, const char ff[])
{
	if (condition)
	{
		cerr << endl << " *=*=*=*=*=*=* ERROR *=*=*=*=*=*=* " << endl<<endl;
		cerr<<"In function: "<<string(ff)<<endl<<endl;
		cerr << error_msg <<endl;
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* " << endl;
		exit(error_code);
	}
}

//
//  dcMatrix.h
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef dcMatrix_h
#define dcMatrix_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <cstdlib>


using namespace std;

class dcMatrix
{
public:
	unsigned long nbRows;	//Dimensions
	unsigned long nbCols;
    
	vector<double> val;	//Value for each element
    
	double & operator()(unsigned long,unsigned long);
	double & operator[](unsigned long);
    
    unsigned long getNbRows()     const {return nbRows;}
    unsigned long getNbCols()     const {return nbCols;}
    
    
    double  get_val(unsigned long i,unsigned long j) const {return val[nbCols*i+j];}
    
    
    // Constructors
	
	//dcMatrix(){nbRows=0;nbCols=0;}
	dcMatrix(){}
	
	dcMatrix(unsigned long h,unsigned long l) 
	{
		nbRows=h; nbCols=l; 
		vector<double> tmp(l*h,0.0);
		val=tmp;
	}
	
	dcMatrix(unsigned long n) 
	{
		nbRows=n; nbCols=n;
		vector<double> tmp(n*n, 0.0);
		val=tmp;
	}
	
    dcMatrix(string pathFile);
	
	dcMatrix(vector<double> v);	// dcMatrix (n,1) from a single vector(n)
	
	
    
    void    resize(unsigned long n)  
	{
		nbRows=n;nbCols=n; 
		val.resize(n*n);
	}
    
	void    resize(unsigned long nRows, unsigned long nCols)  
	{
		nbRows = nRows;
        nbCols = nCols;
		this->val.resize(nRows*nCols);
	}
    
	void    display() const;
        
	void    RandomInit();
    
    // Files operations
	void    FromFile(const char*);
    void    FromFile(string);
	void	FromFile_Rows(string fileName, unsigned long nrow);
	
    void    WriteToFile(string);
	void    WriteToFileCSV(string);
	void    WriteToFileCSV(string fileName, vector<string> header);
	
	
    // Operations on embedded vectors
    vector<double>  extractColumn(unsigned long j_col);
	vector<double>	extractRow(unsigned long i_row);
	
	void            addRowVector(const vector<double>& v);
	void            addRowVector(const vector<unsigned long>& v);
	
    void            addColVector(const vector<double> & v);
    
	void			removeRow(unsigned long i_row);	// removes row 'i_row' and resize dcMatrix
	void			removeCol(unsigned long j_col);	// removes column 'j_col' and resize dcMatrix

	
	
	// Extract the row #i of the matrix
	// "i" is calculated such that it is
	// the smallest element of column "j_col"
  	vector<double>	extractRow_cond_minElement(unsigned long j_col);
	
    // Operations on elements
    
    void    setAllValues(double value);
	void	setValueFromMatrix(dcMatrix M);
	
	double  sumAllElements();
	double  sumLine(unsigned long i);		// sum all elements of line #i
	double  sumColumn(unsigned long j);	// sum all elements of line #i
	
	// conditional sum 
	double	sumColumnIf(unsigned long colToSum, unsigned long colToTest,
						double lowerBound, double upperBound);

	// counts nb of elements which are lower<element<upper  
	unsigned long		countColumnIf(unsigned long colToTest,
						  double lowerBound, double upperBound);
    
    void    setColumnValues(unsigned long colNb, vector<double> v);
	void	setRowValues(unsigned long rowNb_start0, vector<double> v);
	void	setRowValues(unsigned long rowNb_start0, vector<unsigned long> v);
    
    dcMatrix  transpose();
    
    bool    isSymetric();

    double  determinant();
	
	dcMatrix	getMinor(unsigned long row, unsigned long col);
	
	dcMatrix	inverse();
    
    dcMatrix  Cholesky();
    
    double  getMinimumValue();
    double  getMaximumValue();
	
};

dcMatrix operator + (dcMatrix &A,dcMatrix &B);
dcMatrix operator - (dcMatrix &A,dcMatrix &B);
dcMatrix operator * (dcMatrix &A,dcMatrix &B);
dcMatrix operator * (double a,dcMatrix &A);

vector<vector<double> > to_vector_vector(dcMatrix M);

dcMatrix Id(unsigned long n);

dcMatrix power(dcMatrix A,unsigned long n);

dcMatrix cholesky(dcMatrix A);	//renvoi la dcMatrix triangul L tq://si A symetrique,carre //L*transpo(L)=A

double distance_Matrix(dcMatrix A, dcMatrix B, double power);	// Euclidian distance b/w two matrices

dcMatrix rowBind_old(dcMatrix A, dcMatrix B);
//dcMatrix rowBind(dcMatrix A, dcMatrix B);
dcMatrix rowBind(const dcMatrix & A, const dcMatrix & B);

dcMatrix	transpo(dcMatrix A);        // FIX ME : to delete if not used elsewhere
double	Det(dcMatrix A);           // FIX ME : to delete if not used elsewhere
unsigned long		test_sym(dcMatrix A);       // FIX ME : to delete if not used elsewhere



#endif

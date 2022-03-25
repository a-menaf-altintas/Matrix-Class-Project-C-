// This is the header file for the dense real matrix class

#ifndef  _MATRIX_REAL_HPP_
#define  _MATRIX_REAL_HPP_
#include <iostream>
#include <string>
#include <iomanip>
#include <omp.h> // Parallel for loop by using openMP for matrix maultiplication
// STL Containers of the Vector and Array can be used. I will use std::vector
#include <vector> 
//#include <valarray>



template <class T>  // will take care of data structures
class RMatrix{

    private:
        unsigned nrows; // number of row index can only have positive values
        unsigned ncols; // number of coulmn index can only have positive values
        std::vector<std::vector<T> > matrix; // Should be a space between arrows > > Do not confuse with insertion operator
    public:
    //Constructors
        RMatrix(unsigned numberOfRows, unsigned numberOfColumns, const T &value); // Default constructor
        RMatrix(const RMatrix<T> &rhs);  // Copy constructor 
        virtual ~RMatrix();


    // Getters for the number of the row and column                                                                                                                                                                                            
        unsigned get_nrows() const;
        unsigned get_ncols() const;

    // Print the matrix
    void print_matrix();


    // Access the to the  elements  of the matrix                                                                                                                                                                                             
    T &operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns); // Can be modified
    const T &operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns) const; // Read only   


    // Operator overloading
    RMatrix<T> &operator=( const RMatrix<T> &rhs); // matrix_lhs = matrix_rhs

    // Basic binary aritmethic operation of a MATRIX with another MATRIX 
    // Operations Overloading
    template <class Y> friend RMatrix<Y> operator+(const RMatrix<Y> &lhs, const RMatrix<Y> &rhs); //matrix3 = matrix1 + matrix2
    template <class Y> friend RMatrix<Y> &operator+=( RMatrix<Y> &lhs, const RMatrix<Y> &rhs); //matrix1 = matrix1 + matrix2
    
    template <class Y> friend RMatrix<Y> operator-(const RMatrix<Y> &lhs, const RMatrix<Y> &rhs); //matrix3 = matrix1 - matrix2
    template <class Y> friend RMatrix<Y> &operator-=( RMatrix<Y> &lhs, const RMatrix<Y> &rhs); //matrix1 = matrix1 - matrix2
    
    template <class Y> friend RMatrix<Y> operator*(const RMatrix<Y> &lhs, const RMatrix<Y> &rhs); //matrix3 = matrix1 * matrix2
    template <class Y> friend RMatrix<Y> &operator*=( RMatrix<Y> &lhs, const RMatrix<Y> &rhs); //matrix1 = matrix1 * matrix2

    

    // Basic binary aritmethic operation of a MATRIX with a scalar 
    // Operations Overloading

    RMatrix<T> operator+(const T &rhs); // matrix2 = matrix1 + scalar
    template <class Y> friend RMatrix<Y> operator+(const Y &lhs, const RMatrix<Y> &rhs); // matrix2 = scalar + matrix1 

    RMatrix<T> operator-(const T &rhs); // matrix2 = matrix1 - scalar 
    template <class Y> friend RMatrix<Y> operator-(const Y &lhs, const RMatrix<Y> &rhs); //matrix2 =  scalar - matrix1
    
    RMatrix<T> operator*(const T &rhs); // matrix2 = matrix1 * scalar
    template <class Y> friend RMatrix<Y> operator*(const Y &lhs, const RMatrix<Y> &rhs); //matrix2 =  scalar * matrix1


    RMatrix<T> operator/(const T &rhs); // matrix2 = matrix1 / scalar 

    // Basic binary aritmethic operation of a MATRIX with a vector
    // Vector can be only multiplied from right
    // Operations Overloading
    template <class Y> friend std::vector<Y> operator*(const RMatrix<Y> &lhs, const std::vector<Y> &rhs); // matrix2 = matrix1 * vector
};


/************************< IMPLEMENTATION OF CLASS ATRIBUTES FUCNTIONS AND OPERATORS >***********************************************/

// Implementation of default constructor
template <class T>  // will take care of data structures                                                                                                                                                      
RMatrix<T>::RMatrix(unsigned numberOfRows, unsigned numberOfColumns, const T &element) {
  matrix.resize(numberOfRows);
  for (unsigned i=0; i<matrix.size(); i++) {
    matrix[i].resize(numberOfColumns, element);
  }
  nrows = numberOfRows;
  ncols = numberOfColumns;
}

// Copy Constructor                                                                                                                                                           
template<class T>
RMatrix<T>::RMatrix(const RMatrix<T> &rhs) {
  matrix = rhs.matrix;
  nrows = rhs.get_nrows();
  ncols = rhs.get_ncols();
}

// Destructor                                                                                                                                                       
template<class T>
RMatrix<T>::~RMatrix() {}



// Get the number of rows of the matrix                                                                                                                                       
template<class T>
unsigned RMatrix<T>::get_nrows() const {
  return this->nrows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<class T>
unsigned RMatrix<T>::get_ncols() const {
  return this->ncols;
}

// A function to print the matrix.
template <class T>
void RMatrix<T>::print_matrix()
{
	unsigned nRows = this->get_nrows();
	unsigned nCols = this->get_ncols();

  std::cout << std::endl;
	for (unsigned row = 0; row<nRows; row++)
  {
	  for (unsigned col = 0; col<nCols; col++)
    {
	    std::cout <<   this->matrix[row][col] << "  ";
    }
	std::cout << std::endl;
	}
    
}




// Access the to the  elements  of the matrix 
template<class T>
T &RMatrix<T>::operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns) {
  if (numberOfRows >= this->nrows || numberOfColumns >= this->ncols){
    std::cerr<<"\n\nOut of range error!\n\n" <<
    "Index is out of range!\n\n"<<std::endl;
    exit(1);
  }
  return this->matrix[numberOfRows][numberOfColumns];
}

template<class T>
const T &RMatrix<T>::operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns) const {
  if (numberOfRows >= this->nrows || numberOfColumns >= this->ncols){
    std::cerr<<"\n\nOut of range error!\n\n" <<
    "Index is out of range!\n\n"<<std::endl;
    exit(1);
  }
  return this->matrix[numberOfRows][numberOfColumns];
}


// Assignment Operator                                                                                                                                                        
template<class T>
RMatrix<T> &RMatrix<T>::operator=(const RMatrix<T> &rhs) {
  if (this == &rhs)
    return *this;

  unsigned newNumberOfRows = rhs.get_nrows();
  unsigned newNumberOfColumns = rhs.get_ncols();

  matrix.resize(newNumberOfRows);
  for (unsigned i=0; i<matrix.size(); i++) {
    matrix[i].resize(newNumberOfColumns);
  }

  for (unsigned i=0; i<newNumberOfRows; i++) {
    for (unsigned j=0; j<newNumberOfColumns; j++) {
      matrix[i][j] = rhs(i, j);
    }
  }
  nrows = newNumberOfRows;
  ncols = newNumberOfColumns;

  return *this;
}

// Addition operation of a MATRIX with another MATRIX: matrix3 = matrix1 + matrix2                                                                                                                                                 
template<class T>
RMatrix<T> operator+(const RMatrix<T> &lhs, const RMatrix<T> &rhs){
  if (rhs.get_nrows() != lhs.get_nrows() || rhs.get_ncols() != lhs.get_ncols()){
    std::cerr<<"\n\nDimension error while addition of two matrices! << ADDITION: C = A + B >> \n\n" <<
    "Dimensions of A matrix are not equal to dimensions of B matrix!\n\n"<<std::endl;
    exit(1);}
  
  unsigned nrows = rhs.get_nrows();
  unsigned ncols = rhs.get_ncols();
  RMatrix<T> new_matrix(nrows, ncols, 0.0); // A new matrix will be created after addition

  for (unsigned i=0; i<nrows; i++) { 
    for (unsigned j=0; j<ncols; j++) {
      new_matrix(i,j) = lhs(i,j) + rhs(i,j);
    }
  }

  return new_matrix;
}





// Addition operation of a MATRIX with another MATRIX: matrix1 = matrix1 + matrix2: matrix1 += matrix2
template<class T>
RMatrix<T> &operator+=( RMatrix<T> &lhs, const RMatrix<T> &rhs){
  if (rhs.get_nrows() != lhs.get_nrows()  ||  rhs.get_ncols() != lhs.get_ncols()){
    std::cerr<<"\n\nDimension error while addition of two matrices!  << ADDITION: A +=  B >> \n\n" <<
    "Dimensions of A matrix are not equal to dimensions of B matrix!\n\n"<<std::endl;
    exit(1);}
    unsigned nrows = rhs.get_nrows();
    unsigned ncols = rhs.get_ncols();

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        lhs(i,j) = lhs(i,j) + rhs(i,j);
    }
  }

  return lhs;
}


// Subtraction operation of a MATRIX with another MATRIX: matrix3 = matrix1 - matrix2                                                                                                                                                 
template<class T>
RMatrix<T> operator-(const RMatrix<T> &lhs, const RMatrix<T> &rhs){
  if (rhs.get_nrows() != lhs.get_nrows()  ||  rhs.get_ncols() != lhs.get_ncols()){
    std::cerr<<"\n\nDimension error while subtraction of two matrices! << SUBTRACTION: C = A - B >> \n\n" <<
    "Dimensions of A matrix are not equal to dimensions of B matrix!\n\n"<<std::endl;
    exit(1);}
  
  unsigned nrows = rhs.get_nrows();
  unsigned ncols = rhs.get_ncols();
  RMatrix<T> new_matrix(nrows, ncols, 0.0); // A new matrix will be created after subtraction

  for (unsigned i=0; i<nrows; i++) { 
    for (unsigned j=0; j<ncols; j++) {
      new_matrix(i,j) = lhs(i,j) - rhs(i,j);
    }
  }

  return new_matrix;
}





// Subtraction operation of a MATRIX with another MATRIX: matrix1 = matrix1 - matrix2: matrix1 -= matrix2
template<class T>
RMatrix<T> &operator-=( RMatrix<T> &lhs, const RMatrix<T> &rhs){
  if (rhs.get_nrows() != lhs.get_nrows()  ||  rhs.get_ncols() != lhs.get_ncols()){
    std::cerr<<"\n\nDimension error while subtraction of two matrices!  << SUBTRACTION: A -=  B >> \n\n" <<
    "Dimensions of A matrix are not equal to dimensions of B matrix!\n\n"<<std::endl;
    exit(1);}
    unsigned nrows = rhs.get_nrows();
    unsigned ncols = rhs.get_ncols();

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        lhs(i,j) = lhs(i,j) - rhs(i,j);
    }
  }

  return lhs;
}


// Multiplication operation of a MATRIX with another MATRIX: matrix3 = matrix1 * matrix2                                                                                                                                                 
template<class T>
RMatrix<T> operator*(const RMatrix<T> &lhs, const RMatrix<T> &rhs) {
  if (rhs.get_nrows() != lhs.get_ncols()){
    std::cerr<<"\n\nDimension error while multiplication of two matrices! << MULTIPLICATION: C = A * B >> \n\n" <<
    "Number of colmns of A matrix must be equal to number of rows of B matrix!\n\n"<<std::endl;
    exit(1);}
    unsigned ncols = rhs.get_ncols();
    unsigned nrows = lhs.get_nrows();
    unsigned nrows_rhs = rhs.get_nrows();

    RMatrix<T> new_matrix(nrows, ncols, 0.0); // A new matrix will be created after multiplication

 
    #pragma omp parallel for // Outer loop is parallelized by using all cores in the CPU
    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) 
        for (unsigned k=0; k<nrows_rhs; k++) {
        new_matrix(i,j) = new_matrix(i,j) + lhs(i,k) * rhs(k,j);
    }
  }

  return new_matrix;
}

// Multiplication operation of a MATRIX with another MATRIX: matrix1 = matrix1 * matrix2: matrix1 *= matrix2
template<class T>
RMatrix<T> &operator*= ( RMatrix<T> &lhs, const RMatrix<T> &rhs) {
  if ( (rhs.get_nrows() != lhs.get_ncols() ) ||  ( rhs.get_ncols() != rhs.get_nrows() )){
    std::cerr<<"\n\nDimension error while multiplication of two matrices! << MULTIPLICATION: A *= B >> \n\n" <<
    "Number of colmns of A matrix must be equal to number of rows of B matrix!\n\n"<<
    "Matrix B should be a square matrix!\n\n"<<std::endl;
    exit(1);}

    lhs = lhs * rhs; // This operation is already defined before
    return lhs;
}



// Addition operation of a MATRIX with a scalar: matrix1 = matrix1 + scalar
template<class T>
RMatrix<T> RMatrix<T>::operator+(const T &rhs) {
    RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after addition

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = this->matrix[i][j] + rhs;
    }
  }

  return new_matrix;
}


// Addition operation of a scalar with a MATRIX: matrix1 = scalar + matrix1
template<class T>
RMatrix<T> operator+(const T &lhs, const RMatrix<T> &rhs) {
    unsigned nrows = rhs.nrows;
    unsigned ncols = rhs.ncols;
    RMatrix<T> new_matrix(nrows, ncols, 0.0); // A new matrix will be created after addition

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = lhs + rhs(i,j);
    }
  }

  return new_matrix;
}

// Subtraction operation of a MATRIX with a scalar: matrix1 = matrix1 - scalar
template<class T>
RMatrix<T> RMatrix<T>::operator-(const T &rhs) {
    RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after subtraction

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = this->matrix[i][j] - rhs;
    }
  }

  return new_matrix;
}

// Subtraction operation of a scalar  with a MATRIX : matrix1 =  scalar - matrix1
template<class T>
RMatrix<T> operator-(const T &lhs, const RMatrix<T> &rhs) {
    unsigned nrows = rhs.nrows;
    unsigned ncols = rhs.ncols;
    RMatrix<T> new_matrix(nrows, ncols, 0.0); // A new matrix will be created after addition

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = lhs - rhs(i,j);
    }
  }

  return new_matrix;
}




// Multiplication operation of a MATRIX with a scalar: matrix1 = matrix1 * scalar
template<class T>
RMatrix<T> RMatrix<T>::operator*(const T &rhs) {
    RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after subtraction

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = this->matrix[i][j] * rhs;
    }
  }

  return new_matrix;
}

// Multiplication operation of a scalar  with a MATRIX : matrix1 =  scalar * matrix1
template<class T>
RMatrix<T> operator*(const T &lhs, const RMatrix<T> &rhs) {
    unsigned nrows = rhs.nrows;
    unsigned ncols = rhs.ncols;
    RMatrix<T> new_matrix(nrows, ncols, 0.0); // A new matrix will be created after addition

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = lhs * rhs(i,j);
    }
  }

  return new_matrix;
}


// Division operation of a MATRIX with a scalar: matrix1 = matrix1 / scalar
template<class T>
RMatrix<T> RMatrix<T>::operator/(const T &rhs) {
    RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after subtraction

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_matrix(i,j) = this->matrix[i][j] / (double)rhs;
    }
  }

  return new_matrix;
}

// Multiplication operation of a MATRIX with a vector: matrix1 = matrix1 * vector
// Multiplication will return a vector

template<class T>
std::vector<T> operator*(const RMatrix<T> &lhs, const std::vector<T> &rhs) {
  if (lhs.get_ncols() != rhs.size()){
    std::cerr<<"\n\nDimension error while multiplication of a matrix and vector! << MULTIPLICATION: Vec = Mat *  Vec >> \n\n" <<
    "Number of colmns of a matrix must be equal to number of elements of a vector!\n\n"<<std::endl;
    exit(1);}
  
  unsigned nrows = lhs.get_nrows();
  unsigned ncols = lhs.get_ncols();
  std::vector<T> new_vector(nrows, 0.0); // A new vector will be created after multiplication
  

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        new_vector[i] = lhs(i,j) * rhs[j];
    }
  }

  return new_vector;

}


#endif  // end _REAL_MATRIX_HPP_
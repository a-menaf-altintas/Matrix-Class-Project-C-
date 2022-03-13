// This is the header file for the dense real matrix class

#ifndef  _MATRIX_REAL_HPP_
#define  _MATRIX_REAL_HPP_

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


    // Access the to the  elements  of the matrix                                                                                                                                                                                             
    T &operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns);
    const T &operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns) const;    


    // Operator overloading
    RMatrix<T> &operator=( const RMatrix<T> &rhs); // matrix_lhs = matrix_rhs

    // Basic binary aritmethic operation of a MATRIX with another MATRIX 
    // Operations Overloading
    RMatrix<T> operator+(const RMatrix<T> &rhs); // matrix3 = matrix1 + matrix2
    RMatrix<T> &operator+=(const RMatrix<T> &rhs); // matrix1 = matrix1 + matrix2
    RMatrix<T> operator-(const RMatrix<T> &rhs); // matrix3 = matrix1 - matrix2
    RMatrix<T> &operator-=(const RMatrix<T> &rhs); // matrix1 = matrix1 - matrix2
    RMatrix<T> operator*( const RMatrix<T> &rhs); // matrix3 = matrix1 * matrix2
    RMatrix<T> &operator*=( const RMatrix<T> & rhs); // matrix1 = matrix1 * matrix2
    

    // Basic binary aritmethic operation of a MATRIX with a scalar 
    // Operations Overloading

    RMatrix<T> operator+(const T &rhs); // matrix2 = matrix1 + scalar
    template <class Y> friend RMatrix<Y> operator+(const Y &lhs, const RMatrix<Y> &rhs); // matrix2 = scalar + matrix1 

    RMatrix<T> operator-(const T &rhs); // matrix2 = matrix1 - scalar 
    template <class Y> friend RMatrix<Y> operator-(const Y &lhs, const RMatrix<Y> &rhs); //matrix1 =  scalar - matrix1
    
    RMatrix<T> operator*(const T &rhs); // matrix2 = matrix1 * scalar
    template <class Y> friend RMatrix<Y> operator*(const Y &lhs, const RMatrix<Y> &rhs); //matrix1 =  scalar * matrix1


    RMatrix<T> operator/(const T &rhs); // matrix2 = matrix1 / scalar 

    // Basic binary aritmethic operation of a MATRIX with a vector
    // Vector can be only multiplied from right
    // Operations Overloading
    std::vector<T> operator*(const std::vector<T> &rhs); // matrix2 = matrix1 * vector






    // RMatrix<T> transpose(); // take the transpose of a matrix





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

// (Virtual) Destructor                                                                                                                                                       
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

// Access the to the  elements  of the matrix 
template<class T>
T &RMatrix<T>::operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns) {
  return this->matrix[numberOfRows][numberOfColumns];
}

template<class T>
const T &RMatrix<T>::operator()(const unsigned &numberOfRows, const unsigned &numberOfColumns) const {
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
RMatrix<T> RMatrix<T>::operator+(const RMatrix<T> &rhs) {
  RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after addition

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      new_matrix(i,j) = this->matrix[i][j] + rhs(i,j);
    }
  }

  return new_matrix;
}

// Addition operation of a MATRIX with another MATRIX: matrix1 = matrix1 + matrix2: matrix1 += matrix2
template<class T>
RMatrix<T> &RMatrix<T>::operator+=(const RMatrix<T> &rhs) {
    unsigned nrows = rhs.get_nrows();
    unsigned ncols = rhs.get_ncols();

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        this->matrix[i][j] = this->matrix[i][j] + rhs(i,j);
    }
  }

  return *this;
}


// Subtraction operation of a MATRIX with another MATRIX: matrix3 = matrix1 - matrix2                                                                                                                                                 
template<class T>
RMatrix<T> RMatrix<T>::operator-(const RMatrix<T> &rhs) {
  RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after subtraction

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      new_matrix(i,j) = this->matrix[i][j] - rhs(i,j);
    }
  }

  return new_matrix;
}

// Subtraction operation of a MATRIX with another MATRIX: matrix1 = matrix1 - matrix2: matrix1 -= matrix2
template<class T>
RMatrix<T> &RMatrix<T>::operator-= (const RMatrix<T> &rhs) {
    unsigned nrows = rhs.get_nrows();
    unsigned ncols = rhs.get_ncols();

    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) {
        this->matrix[i][j] = this->matrix[i][j] - rhs(i,j);
    }
  }

  return *this;
}


// Multiplication operation of a MATRIX with another MATRIX: matrix3 = matrix1 * matrix2                                                                                                                                                 
template<class T>
RMatrix<T> RMatrix<T>::operator*(const RMatrix<T> &rhs) {
    unsigned nrows = rhs.get_nrows();
    unsigned ncols = rhs.get_ncols();

    RMatrix new_matrix(nrows, ncols, 0.0); // A new matrix will be created after multiplication


    for (unsigned i=0; i<nrows; i++) {
        for (unsigned j=0; j<ncols; j++) 
        for (unsigned k=0; k<nrows; k++) {
        new_matrix(i,j) = new_matrix(i,j) + this->matrix[i][k] * rhs(k,j);
    }
  }

  return new_matrix;
}

// Multiplication operation of a MATRIX with another MATRIX: matrix1 = matrix1 * matrix2: matrix1 *= matrix2
template<class T>
RMatrix<T> &RMatrix<T>::operator*= (const RMatrix<T> &rhs) {
    RMatrix new_matrix = (*this) * rhs; // This operation already defined above
    (*this) = new_matrix;
    return *this;
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






#endif  // end _REAL_MATRIX_HPP_
#if !defined(__Matrix_hdr__)
#define __Matrix_hdr__


class CMatrix
{
public:

    enum MatrixError_t 
    {
        kMatrixNoError,

        kMatrixNotSquare,
        kMatrixDimensionMisfit,

        kMatrixMemAllocError,

        kMatrixIllegalPointer,

        kMatrixIllegalFunctionParam,

        kMatrixUnknownError,

        kMatrixNumErrors
    };
    CMatrix ();
    CMatrix(int iNumRows, int iNumCols);
    virtual ~CMatrix ();
    CMatrix (const CMatrix &other);

    void assign( const CMatrix &other );

    MatrixError_t   init (int iNumRows, int iNumCols);
    MatrixError_t   reset ();

    int             getNumRows () const;
    int             getNumCols () const;

    MatrixError_t   setElement (int iRow, int iCol, float fValue);
    float           getElement (int iRow, int iCol) const;
    MatrixError_t   setRow (int iRow, const float *pfValues, int iNumValues);
    MatrixError_t   getRow (int iRow, float *pfValues, int iNumValues) const;
    MatrixError_t   setCol (int iCol, const float *pfValues, int iNumValues);
    MatrixError_t   getCol (int iCol, float *pfValues, int iNumValues) const;

    MatrixError_t   setZero ();
    MatrixError_t   setRand ();
    MatrixError_t   transpose ();
    MatrixError_t   normalize (int p = 1);

    MatrixError_t   mulByElement(const CMatrix &other);
    MatrixError_t   divByElement(const CMatrix &other);
    MatrixError_t   addByElement(const CMatrix &other);
    MatrixError_t   subByElement(const CMatrix &other);

    CMatrix& operator = (const CMatrix &other);
    CMatrix operator + (const CMatrix &other) const;
    CMatrix operator - (const CMatrix &other) const;
    CMatrix operator * (const CMatrix &other) const;
    bool    operator== (const CMatrix &other) const;
    CMatrix operator + (const float fValue) const;
    CMatrix operator - (const float fValue) const;
    CMatrix operator * (const float fValue) const;

    float getNorm(int p = 1) const;


private:
    const float*    getRow (int iRow) const;
    bool isIndexValid( int iRow, int iCol ) const
    {
        if ((iRow >= m_aiMatrixDimensions[kRow]) || (iCol >= m_aiMatrixDimensions[kCol]) || iRow < 0 || iCol < 0)
            return false;
        return true;
    }

    enum Dimensions_t
    {
        kRow,
        kCol,

        kNumOfDimensions
    };

    int        m_aiMatrixDimensions[kNumOfDimensions];

    float    **m_ppfMatrix;

};

#endif // #if !defined(__Matrix_hdr__)




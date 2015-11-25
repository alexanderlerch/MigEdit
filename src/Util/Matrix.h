#if !defined(__MatrixView_hdr__)
#define __MatrixView_hdr__


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
    virtual ~CMatrix ();
    CMatrix (const CMatrix &other);

    MatrixError_t   initialize (int iNumRows, int iNumColumns);
    MatrixError_t   reset ();

    int             getNumRows () const;
    int             getNumCols () const;

    MatrixError_t   setZero ();

    MatrixError_t   setElement (int iRow, int iCol, float fValue);
    float           getElement (int iRow, int iCol) const;
    MatrixError_t   setRow (int iRow, const float *pfValues, int iNumValues);
    MatrixError_t   getRow (int iRow, float *pfValues, int iNumValues) const;
    MatrixError_t   setCol (int iCol, const float *pfValues, int iNumValues);
    MatrixError_t   getCol (int iCol, float *pfValues, int iNumValues) const;

    CMatrix operator = (const CMatrix &other);
    CMatrix operator + (const CMatrix &other);
    CMatrix operator - (const CMatrix &other);
    CMatrix operator * (const CMatrix &other);
    bool    operator== (const CMatrix &other);
    CMatrix operator + (const float fValue);
    CMatrix operator - (const float fValue);
    CMatrix operator * (const float fValue);


private:
    bool isIndexValid( int iRow, int iCol ) const
    {
        if ((iRow >= m_aiMatrixDimensions[kRow]) || (iCol >= m_aiMatrixDimensions[kCol]) || iRow < 0 || iCol < 0)
            return false;
        return true;
    }
    //void   RowDivC (int iRow, float fValue);
    //void   RowMultC (int iRow, float fValue);
    //void   RowMultC (int iRow, int iColBegin, float fValue);
    //void   RowAdd (int iRow1, int iRow2, float fValue);
    //void   RowSub (int iRow1, int iRow2, float fValue);
    //void   ColDivC (int iCol, float fValue);
    //void   ColMultC (int iCol, float fValue);
    //void   ColMultC (int iCol, int iRowBegin, float fValue);
    //void   ColAdd (int iCol1, int iCol2, float fValue);
    //void   ColAddC (int iCol, float fValue);

    //MatrixError_t HouseholderCornerSubtraction(CMatrix &Smaller);
    //MatrixError_t HouseHolder (CMatrix &Input);
    //MatrixError_t HouseHolderMult (CMatrix &VectorInput);
    //MatrixError_t GetPartialMatrix (CMatrix *pResult, int iStartRow, int iStartCol);
    //MatrixError_t CalcScalarProd (CMatrix &V, float *pfResult);

    //void   SwapRowCol (int iRow, int iCol);
    //MatrixError_t Inverse ();
    //MatrixError_t GetSubMatrix (CMatrix *pResult, int iRow, int iCol);


    enum Dimensions_t
    {
        kRow,
        kCol,

        kNumOfDimensions
    };

    int        m_aiMatrixDimensions[kNumOfDimensions];

    float    **m_ppfMatrix;

};

#endif // #if !defined(__MatrixView_hdr__)




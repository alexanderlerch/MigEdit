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

    enum ActionAppliedTo_t
    {
        kPerRow,
        kPerCol,
        kAll,

        kNumActionAreas
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
    const float*    getRowPtr (int iRow) const;
    MatrixError_t   setCol (int iCol, const float *pfValues, int iNumValues);
    MatrixError_t   getCol (int iCol, float *pfValues, int iNumValues) const;
    CMatrix         getSubMatrix(int iStartRow, int iStartCol, int iEndRow, int iEndCol) const;

    MatrixError_t   setZero ();
    MatrixError_t   setZeroBelowThresh (float fThresh = 0);
    MatrixError_t   setOnes ();
    MatrixError_t   setRand ();
    MatrixError_t   setRowRand (int iRow);
    
    MatrixError_t   transpose_I ();
    MatrixError_t   normalize_I (ActionAppliedTo_t eActionArea = kAll, int p = 1);
    CMatrix         transpose ();
    CMatrix         normalize (ActionAppliedTo_t eActionArea = kAll, int p = 1);

    MatrixError_t   mulByOnes_I(int iNumRows, int iNumCols);
    CMatrix         mulByOnes(int iNumRows, int iNumCols);

    CMatrix&        mulByElement_I(const CMatrix &other);
    CMatrix&        mulC_I(float fC);
    CMatrix&        mulRowC_I(int iRow, float fC);
    CMatrix&        divByElement_I(const CMatrix &other);
    CMatrix&        addByElement_I(const CMatrix &other);
    CMatrix&        addC_I( float fC );
    MatrixError_t   subByElement_I(const CMatrix &other);
    CMatrix         mulByElement(const CMatrix &other) const;
    CMatrix         divByElement(const CMatrix &other) const;
    CMatrix         addByElement(const CMatrix &other) const;
    CMatrix         subByElement(const CMatrix &other) const;

    CMatrix& operator = (const CMatrix &other);
    CMatrix operator + (const CMatrix &other) const;
    CMatrix operator - (const CMatrix &other) const;
    CMatrix operator * (CMatrix &other) const;
    bool    operator== (const CMatrix &other) const;
    CMatrix operator + (const float fValue) const;
    CMatrix operator - (const float fValue) const;
    CMatrix operator * (const float fValue) const;
    static float getKlDivergence( const CMatrix &Mat1, const CMatrix &Mat2 ) 
    {
        CMatrix *pMat  = 0;
        int iNumRows    = Mat1.getNumRows();
        int iNumCols    = Mat1.getNumCols();
        float fResult   = 0;
        float **ppfMat1 = 0;
        float **ppfMat2 = 0;

        assert (iNumRows == Mat2.getNumRows() || iNumCols == Mat2.getNumCols());

        pMat            = (CMatrix*)(&Mat1);
        ppfMat1         = pMat->getMatrixPtr();
        pMat            = (CMatrix*)(&Mat2);
        ppfMat2         = pMat->getMatrixPtr();
        for (int i = 0; i < iNumRows; i++)
        {
            for (int j = 0; j < iNumCols; j++)
            {
                //fResult     += ppfMat1[i][j] * (std::logf(ppfMat1[i][j] + 1e-24F) - std::logf(ppfMat2[i][j] + 1e-24F)) - ppfMat1[i][j] + ppfMat2[i][j];
                fResult     += ppfMat1[i][j] * (std::logf((ppfMat1[i][j] + 1e-24F)/(ppfMat2[i][j] + 1e-24F)) - 1.F) + ppfMat2[i][j];
            }
        }

        return fResult;
    }


    float getNorm(int p = 1) const;
    float getRowNorm(int iRow, int p = 1) const;
    float getColNorm(int iCol, int p = 1) const;
    float getSum(bool bAbs = false) const;
    float getMax(bool bAbs = false) const;
    void dbgPrint2StdOut() const;

private:
    float** getMatrixPtr();
    float getVectorNorm (int iRow = -1, int iCol = -1, int p = 1) const;
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




#include <string>
#include <cassert>

#include "Util.h"
#include "Matrix.h"



CMatrix::CMatrix () :
    m_ppfMatrix(0)
{
    m_aiMatrixDimensions[kRow]  = 0;
    m_aiMatrixDimensions[kCol]  = 0;
}

CMatrix::CMatrix( const CMatrix &other )
{
    initialize(other.getNumRows(), other.getNumCols());
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        other.getRow(i, m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);
}


CMatrix::~CMatrix ()
{
    this->reset ();
}

CMatrix::MatrixError_t CMatrix::reset ()
{
    if (m_ppfMatrix)
    {
        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            if (m_ppfMatrix[i])
                delete [] m_ppfMatrix[i];
            m_ppfMatrix[i]  = 0;
        }
    }
    delete [] m_ppfMatrix;
    m_ppfMatrix = 0;

    m_aiMatrixDimensions[kRow]  = 0;
    m_aiMatrixDimensions[kCol]  = 0;

    return kMatrixNoError;
}


CMatrix::MatrixError_t   CMatrix::initialize (int iNumRows, int iNumCols)
{
    int i;

    assert (iNumRows > 0);
    assert (iNumCols > 0);

    if (iNumRows == m_aiMatrixDimensions[kRow] && iNumCols == m_aiMatrixDimensions[kCol])
    {
        setZero();
        return kMatrixNoError;
    }
    else
        reset();

    m_ppfMatrix = new float* [iNumRows];

    if (!m_ppfMatrix)
        return kMatrixMemAllocError;

    for (i = 0; i < iNumRows; i++)
    {
        m_ppfMatrix[i]  = new float [iNumCols];
        if (!m_ppfMatrix[i])
            return kMatrixMemAllocError;
    }

    m_aiMatrixDimensions[kRow]  = iNumRows;
    m_aiMatrixDimensions[kCol]  = iNumCols;

    setZero();

    return kMatrixNoError;
}

int            CMatrix::getNumRows () const
{
    return m_aiMatrixDimensions[kRow];
}

int            CMatrix::getNumCols () const
{
    return m_aiMatrixDimensions[kCol];
}


CMatrix::MatrixError_t   CMatrix::setElement (int iRow, int iCol, float fValue)
{
    if(!isIndexValid(iRow, iCol))
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);
    assert (m_ppfMatrix[iRow] != 0);
    
    m_ppfMatrix[iRow][iCol] = fValue;

    return kMatrixNoError;
}

float        CMatrix::getElement (int iRow, int iCol) const
{
    if(!isIndexValid(iRow, iCol))
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);
    assert (m_ppfMatrix[iRow] != 0);

    return m_ppfMatrix[iRow][iCol];
}

CMatrix::MatrixError_t   CMatrix::setRow (int iRow, const float *pfValues, int iNumOfValues)
{
    if(!isIndexValid(iRow, 0) || iNumOfValues != m_aiMatrixDimensions[kCol])
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);
    assert (m_ppfMatrix[iRow] != 0);

    CUtil::copyBuff(m_ppfMatrix[iRow], pfValues, iNumOfValues);

    return kMatrixNoError;
}

CMatrix::MatrixError_t   CMatrix::getRow (int iRow, float *pfValues, int iNumOfValues) const
{
    if(!isIndexValid(iRow, 0) || iNumOfValues != m_aiMatrixDimensions[kCol])
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);
    assert (m_ppfMatrix[iRow] != 0);

    CUtil::copyBuff(pfValues, m_ppfMatrix[iRow], iNumOfValues);

    return kMatrixNoError;
}

CMatrix::MatrixError_t   CMatrix::setCol (int iCol, const float *pfValues, int iNumOfValues)
{
    if(!isIndexValid(0, iCol) || iNumOfValues != m_aiMatrixDimensions[kRow])
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);

    for (int i = 0; i < iNumOfValues; i++)
        m_ppfMatrix[i][iCol]    = pfValues[i];

    return kMatrixNoError;
}

CMatrix::MatrixError_t   CMatrix::getCol (int iCol, float *pfValues, int iNumOfValues) const
{
    if(!isIndexValid(0, iCol) || iNumOfValues != m_aiMatrixDimensions[kRow])
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);

    for (int i = 0; i < iNumOfValues; i++)
        pfValues[i] = m_ppfMatrix[i][iCol];

    return kMatrixNoError;
}

CMatrix CMatrix::operator=(const CMatrix &other)
{ 
    CMatrix		MatrixResult = CMatrix(other);

    return MatrixResult;
}

CMatrix CMatrix::operator+(const CMatrix &other)
{ 
    assert(m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols());

    if (m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols())
    { 
        CMatrix		MatrixResult;
        MatrixResult.initialize(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
            {
                MatrixResult.setElement(i, j, m_ppfMatrix[i][j] + other.getElement(i,j));
            }
        }
        return MatrixResult;
    }
    return CMatrix();
}

CMatrix CMatrix::operator-(const CMatrix &other)
{ 
    assert(m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols());

    if (m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols())
    { 
        CMatrix		MatrixResult;
        MatrixResult.initialize(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
            {
                MatrixResult.setElement(i, j, m_ppfMatrix[i][j] - other.getElement(i,j));
            }
        }
        return MatrixResult;
    }
    return CMatrix();
}

CMatrix CMatrix::operator*(const CMatrix &other)
{ 
    assert(m_aiMatrixDimensions[kCol] == other.getNumRows());

    if (m_aiMatrixDimensions[kCol] != other.getNumRows())
    { 
        CMatrix		MatrixResult;
        MatrixResult.initialize(m_aiMatrixDimensions[kRow], other.getNumCols());

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            for (int j = 0; j < other.getNumCols (); j++)
            {
                float fResult = 0;
                for (int k = 0; k < m_aiMatrixDimensions[kRow]; k++)
                    fResult     += m_ppfMatrix[i][k] * other.getElement(k,j);
                MatrixResult.setElement(i, j, fResult);
            }
        }
        return MatrixResult;
    }
    return CMatrix();
}

bool   CMatrix::operator==(const CMatrix &other)
{
    float    *pfTempValues = 0;
    if (&other == this)
        return true;
    if (m_aiMatrixDimensions[kRow] != other.getNumRows () || m_aiMatrixDimensions[kCol] != other.getNumCols ())
        return false;

    pfTempValues    = new float [m_aiMatrixDimensions[kCol]];

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
    {
        MatrixError_t   rErr;
        rErr    = other.getRow (i, pfTempValues, m_aiMatrixDimensions[kCol]);
        if (rErr)
        {
            delete [] pfTempValues;
            return false;
        }
        if (!CUtil::isBuffEqual(pfTempValues, m_ppfMatrix[i], m_aiMatrixDimensions[kCol]))
        {
            delete [] pfTempValues;
            return false;
        }
    }
    delete [] pfTempValues;
    return true;
}

CMatrix CMatrix::operator+(const float fValue)
{ 
    CMatrix		MatrixResult;
    MatrixResult.initialize(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
    {
        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
        {
            MatrixResult.setElement(i, j, m_ppfMatrix[i][j] * fValue);
        }
    }
    return MatrixResult;
}

//CMatrix::MatrixError_t   CMatrix::CalcMultiply (CMatrix *pSrc)
//{
//    float    fValue;
//    int        iNumNewColumns  = pSrc->GetNumColumns (),
//        iNumOldColumns  = pSrc->GetNumRows (),
//        iNumRows        = this->GetNumRows ();
//
//    CMatrix     MatrixResult;
//    CMatrix::MatrixError_t rErr;
//
//    if (this->GetNumColumns () != iNumOldColumns)
//        return kMatrixDimensionMisfit;
//
//    rErr    = MatrixResult.SetDimensions (iNumRows, iNumNewColumns);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    // e.g.
//    // [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
//    // [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
//    //             [K][L]
//    //
//    for (int i = 0; i < iNumRows; i++)
//    {
//        for (int j = 0; j < iNumNewColumns; j++)
//        {
//            fValue  = 0.0;
//            for (int k = 0; k < iNumOldColumns ; k++)
//                fValue += this->GetMatrixElement(i, k) * pSrc->GetMatrixElement(k, j) ;
//
//            MatrixResult.SetMatrixElement(i, j, fValue) ;
//        }
//    }
//
//    return this->Copy (MatrixResult);
//}
//
//CMatrix::MatrixError_t CMatrix::GetMultiply (CMatrix *pDst, CMatrix *pRight)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    rErr    = pDst->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pDst->CalcMultiply (pRight);
//}
//
//CMatrix::MatrixError_t   CMatrix::CalcAdd (CMatrix *pSrc)
//{
//    // now add in the other matrix
//    int    iNumRows = pSrc->GetNumRows (),
//        iNumCols = pSrc->GetNumColumns ();
//    // first check for a valid addition operation
//    if ((this->GetNumColumns () != iNumCols) || (this->GetNumRows () != iNumRows))
//        return kMatrixDimensionMisfit;
//
//    for (int i = 0 ; i < iNumRows ; i++)
//    {
//        for (int j = 0 ; j < iNumCols ; j++)
//            this->Add2MatrixElement (i, j, pSrc->GetMatrixElement(i, j));
//    }
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::GetAdd (CMatrix *pDst, CMatrix *pRight)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    if (m_aiMatrixDimensions[kRow] != m_aiMatrixDimensions[kCol])
//        return kMatrixNotSquare;
//
//    rErr    = pDst->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pDst->CalcAdd (pRight);
//}
//
//CMatrix::MatrixError_t   CMatrix::CalcSub (CMatrix *pSrc)
//{
//    // now add in the other matrix
//    int    iNumRows = pSrc->GetNumRows (),
//        iNumCols = pSrc->GetNumColumns ();
//
//    // first check for a valid addition operation
//    if ((this->GetNumColumns () != iNumCols) || (this->GetNumRows () != iNumRows))
//        return kMatrixDimensionMisfit
//        ;
//
//    for (int i = 0 ; i < iNumRows ; i++)
//    {
//        for (int j = 0 ; j < iNumCols ; j++)
//            this->Add2MatrixElement (i, j, - pSrc->GetMatrixElement(i, j));
//    }
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::GetSub (CMatrix *pDst, CMatrix *pRight)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    if (m_aiMatrixDimensions[kRow] != m_aiMatrixDimensions[kCol])
//        return kMatrixNotSquare;
//
//    rErr    = pDst->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pDst->CalcSub (pRight);
//}
//
//CMatrix::MatrixError_t   CMatrix::CalcDiagFromVec ()
//{
//    int        iNumCols    = m_aiMatrixDimensions[kCol],
//        iNumRows    = m_aiMatrixDimensions[kRow],
//        iDimensions;
//    bool       bRowVec     = (iNumRows > iNumCols) ? true : false;
//
//    CMatrix     MatrixResult;
//    CMatrix::MatrixError_t rErr;
//
//
//    if (iNumRows != 1 && iNumCols != 1)
//        return kMatrixDimensionMisfit;
//
//    iDimensions = std::max (iNumRows, iNumCols);
//    rErr    = MatrixResult.SetDimensions (iDimensions, iDimensions);
//    if (rErr != kMatrixNoError)
//        return rErr;
//    rErr    = MatrixResult.SetMatrixZero ();
//
//    for (int i = 0; i < iDimensions; i++)
//    {
//        int    iRowIdx = (bRowVec) ? i : 0,
//            iColIdx = (bRowVec) ? 0 : i;
//        rErr    = MatrixResult.SetMatrixElement (i, i, m_ppfMatrix[iRowIdx][iColIdx]);
//        if (rErr != kMatrixNoError)
//            return rErr;
//    }
//
//    return this->Copy (MatrixResult);
//}
//CMatrix::MatrixError_t   CMatrix::CalcVecFromDiag ()
//{
//    int        iNumCols    = m_aiMatrixDimensions[kCol],
//        iNumRows    = m_aiMatrixDimensions[kRow],
//        iDimensions;
//    //bool       bRowVec     = (iNumRows > iNumCols) ? true : false;
//
//    CMatrix     MatrixResult;
//    CMatrix::MatrixError_t rErr;
//
//
//    iDimensions = std::min (iNumRows, iNumCols);
//    rErr    = MatrixResult.SetDimensions (iDimensions, 1);
//    if (rErr != kMatrixNoError)
//        return rErr;
//    rErr    = MatrixResult.SetMatrixZero ();
//
//    for (int i = 0; i < iDimensions; i++)
//    {
//        rErr    = MatrixResult.SetMatrixElement (i, 0, m_ppfMatrix[i][i]);
//        if (rErr != kMatrixNoError)
//            return rErr;
//    }
//
//    return this->Copy (MatrixResult);
//}
//
//CMatrix::MatrixError_t CMatrix::GetDiagFromVec (CMatrix *pResult)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    rErr    = pResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pResult->CalcDiagFromVec ();
//}
//
//CMatrix::MatrixError_t CMatrix::GetVecFromDiag (CMatrix *pResult)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    rErr    = pResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pResult->CalcVecFromDiag ();
//}
//
//CMatrix::MatrixError_t CMatrix::CalcTriu ()
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//    {
//        int iStop = std::min (i, m_aiMatrixDimensions[kCol]);
//        for (int j = 0; j < iStop; j++)
//            m_ppfMatrix[i][j]   = 0;
//    }
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::GetTriu (CMatrix *pResult)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    rErr    = pResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pResult->CalcTriu ();
//}
//
//CMatrix::MatrixError_t CMatrix::CalcTril ()
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//    {
//        for (int j = i+1; j < m_aiMatrixDimensions[kCol]; j++)
//            m_ppfMatrix[i][j]   = 0;
//    }
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::GetTril (CMatrix *pResult)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    rErr    = pResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pResult->CalcTril ();
//}
//
//CMatrix::MatrixError_t   CMatrix::CalcCovariant ()
//{
//    CMatrix         Tmp;
//    CMatrix::MatrixError_t   rErr;
//
//
//    // remove mean
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//    {
//        float dMean  = this->CalcColMean (j);
//        this->ColAddC (j, -dMean);
//    }
//
//    rErr    = Tmp.Copy (*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    rErr    = this->CalcTransposition ();
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    rErr    = this->CalcMultiply (&Tmp);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    // normalize
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        this->RowMultC (i, 1.0/(m_aiMatrixDimensions[kRow]-1));
//
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::GetCovariate (CMatrix *pResult)
//{
//    CMatrix::MatrixError_t   rErr;
//
//    if (m_aiMatrixDimensions[kRow] != m_aiMatrixDimensions[kCol])
//        return kMatrixNotSquare;
//
//    rErr    = pResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pResult->CalcCovariant ();
//}
//
//// a = uwv, u replaces a (*this)
//CMatrix::MatrixError_t CMatrix::CalcSVD (CMatrix *pWResult, CMatrix *pVResult, int iMaxIterations)
//{
//    //    CMatrix::MatrixError_t rErr;
//    //    int        i,j,k,l;
//    //
//    //    float    *pdRV1; 
//    //    float    dTmp, 
//    //                dF,
//    //                dH, 
//    //                dS,
//    //                dNorm       = 0.0,
//    //                dG          = 0.0,
//    //                dScale      = 0.0;
//    //    int        iNumRows   = m_aiMatrixDimensions[kRow],
//    //                iNumCols   = m_aiMatrixDimensions[kCol];
//    //
//    //    float **a; 
//    //
//    //    if (iNumRows < iNumCols) 
//    //        return kMatrixNoSVD;
//    //
//    //
//    //    a= new float* [iNumRows];
//    //    for (i = 0; i < iNumRows; i++)
//    //    {
//    //        float aCopy[10];
//    //        a[i]= new float [iNumCols];
//    //        this->GetMatrixRow (i, aCopy, iNumCols);
//    //        memcpy (&a[i][0], aCopy,sizeof(float)*iNumCols);
//    //    }
//    //
//    //
//    //    pdRV1       = new float [iNumCols]; 
//    //    pWResult->SetDimensions (1, iNumCols);
//    //    pVResult->SetDimensions (iNumRows, iNumCols);
//    //
//    //
//    //    for (i=0;i<iNumCols;i++) 
//    //    {
//    //        l=i+1;
//    //        pdRV1[i]=dScale*dG;
//    //        dG=dS=dScale=0.0;
//    //        if (i < iNumRows) 
//    //        {
//    //            for (k=i;k<iNumRows;k++) 
//    //                dScale += fabs (this->GetMatrixElement (k,i));
//    ////                dScale += fabs(a[k][i]);
//    //            if (dScale) 
//    //            {
//    //                this->ColMultC (i, i, 1/dScale);
//    //                for (k=i;k<iNumRows;k++) 
//    //                {
//    //                    a[k][i] /= dScale;
//    //                    dTmp    = this->GetMatrixElement (k, i);
//    ////                    dS     += dTmp * dTmp;
//    //                    dS += a[k][i]*a[k][i];
//    //                }
//    ////                dF  = this->GetMatrixElement (i, i);
//    //                dF=a[i][i];
//    //                dG = -MHLP_SIGN(sqrt(dS),dF);
//    //                dH=dF*dG-dS;
//    //                this->SetMatrixElement (i, i, dF-dG);
//    //                a[i][i]=dF-dG;
//    //                if (i != iNumCols) 
//    //                {
//    //                    for (j=l;j<iNumCols;j++) 
//    //                    {
//    //                        for (dS=0.0,k=i;k<iNumRows;k++) 
//    //                            dS += a[k][i]*a[k][j];
//    //                        dF=dS/dH;
//    //                        for (k=i;k<iNumRows;k++) 
//    //                        {
//    //                            this->Add2MatrixElement (k, j, dF * a[k][i]);
//    //                            a[k][j] += dF*a[k][i];
//    //                        }
//    //                    }
//    //                }
//    //                ColMultC (i, i, dScale);
//    //                for (k=i;k<iNumRows;k++) 
//    //                    a[k][i] *= dScale;
//    //            }
//    //        }
//    //        pWResult->SetMatrixElement (0, i, dScale*dG);
//    //        dG=dS=dScale=0.0;
//    //        if (i <= iNumRows && i != iNumCols) 
//    //        {
//    //            for (k=l;k<iNumCols;k++) 
//    //                dScale += fabs(a[i][k]);
//    //            if (dScale) 
//    //            {
//    //                RowMultC (i, l, 1/dScale);
//    //                for (k=l;k<iNumCols;k++) 
//    //                {
//    //                    a[i][k] /= dScale;
//    //                    dS += a[i][k]*a[i][k];
//    //                }
//    //                dF=a[i][l];
//    //                dG = -MHLP_SIGN(sqrt(dS),dF);
//    //                dH=dF*dG-dS;
//    //                this->SetMatrixElement (i, l, dF-dG);
//    //                a[i][l]=dF-dG;
//    //                for (k=l;k<iNumCols;k++) 
//    //                    pdRV1[k]=a[i][k]/dH;
//    //                if (i != iNumRows) 
//    //                {
//    //                    for (j=l;j<iNumRows;j++) 
//    //                    {
//    //                        for (dS=0.0,k=l;k<iNumCols;k++) 
//    //                            dS += a[j][k]*a[i][k];
//    //                        for (k=l;k<iNumCols;k++) 
//    //                        {
//    //                            this->Add2MatrixElement (j, k, dS*pdRV1[k]);
//    //                            a[j][k] += dS*pdRV1[k];
//    //                        }
//    //                    }
//    //                }
//    //                RowMultC (i, l, dScale);
//    //                for (k=l;k<iNumCols;k++) 
//    //                    a[i][k] *= dScale;
//    //            }
//    //        }
//    //        {
//    //            /* MAX may not be efficient about the number of times
//    //            it evaluates the fabs sum. This will fix it. */
//    //            float arg2;
//    //            arg2= fabs (pWResult->GetMatrixElement (0, i)) + fabs (pdRV1[i]);
//    //            dNorm=std::max(dNorm, arg2);
//    //        }
//    //    }
//    //    for (i=iNumCols-1;i>=0;i--) 
//    //    {
//    //        if (i < iNumCols) 
//    //        {
//    //            if (dG) 
//    //            {
//    //                for (j=l;j<iNumCols;j++)
//    //                {
//    //                    dTmp    = a[i][j]/a[i][l];
//    //                    pVResult->SetMatrixElement (j, i, dTmp / dG);
//    //                }
//    //                for (j=l;j<iNumCols;j++) 
//    //                {
//    //                    for (dS=0.0,k=l;k<iNumCols;k++) 
//    //                        dS += a[i][k]* pVResult->GetMatrixElement (k, j);
//    //                    for (k=l;k<iNumCols;k++) 
//    //                        pVResult->Add2MatrixElement (k, j, dS * pVResult->GetMatrixElement (k, i));
//    //                }
//    //            }
//    //            for (j=l;j<iNumCols;j++) 
//    //            {
//    //                pVResult->SetMatrixElement (i, j, 0);
//    //                pVResult->SetMatrixElement (j, i, 0);
//    //            }
//    //        }
//    //        pVResult->SetMatrixElement (i, i, 1);
//    //        dG= pdRV1[i];
//    //        l=i;
//    //    }
//    //    for (i=iNumCols-1;i>=0;i--) 
//    //    {
//    //        l=i+1;
//    //        dG= pWResult->GetMatrixElement (0, i);
//    //        if (i < iNumCols)
//    //        {
//    //            RowMultC (i, l, 0);
//    //            for (j=l;j<iNumCols;j++) 
//    //                a[i][j]=0.0;
//    //        }
//    //        if (dG) 
//    //        {
//    //            dG=1.0/dG;
//    //            if (i != iNumCols) 
//    //            {
//    //                for (j=l;j<iNumCols;j++) 
//    //                {
//    //                    for (dS=0.0,k=l;k<iNumRows;k++) 
//    //                        dS += a[k][i]*a[k][j];
//    //                    dF=(dS/a[i][i])*dG;
//    //                    for (k=i;k<iNumRows;k++) 
//    //                    {
//    //                        this->Add2MatrixElement (k, j, dF * a[k][i]);
//    //                        a[k][j] += dF*a[k][i];
//    //                    }
//    //                }
//    //            }
//    //            ColMultC (i, i, dG);
//    //            for (j=i;j<iNumRows;j++) 
//    //                a[j][i] *= dG;
//    //        } 
//    //        else 
//    //        {
//    //            ColMultC (i, i, 0);
//    //            for (j=i;j<iNumRows;j++) 
//    //                a[j][i]=0.0;
//    //        }
//    //        this->Add2MatrixElement (i, i, 1);
//    //        ++a[i][i];
//    //    }
//    //    for (k=iNumCols-1;k>=0;k--) 
//    //    {
//    //        for (int iIter=0;iIter<iMaxIterations;iIter++) 
//    //        {
//    //            float    dC, dX, dY, dZ;
//    //            int        iNM     = 0,
//    //                        iFlag   = 1;
//    //
//    //            for (l=k;l>=0;l--) 
//    //            {
//    //                iNM=l-1;
//    //                if (fabs(pdRV1[l])+dNorm == dNorm) 
//    //                {
//    //                    iFlag=0;
//    //                    break;
//    //                }
//    //                if (fabs (pWResult->GetMatrixElement (0, iNM)) + dNorm == dNorm) 
//    //                    break;
//    //            }
//    //            if (iFlag) 
//    //            {
//    //                dC=0.0;
//    //                dS=1.0;
//    //                for (i=l;i<=k;i++) 
//    //                {
//    //                    dF=dS*pdRV1[i];
//    //                    if (fabs(dF)+dNorm != dNorm) 
//    //                    {
//    //                        dG= pWResult->GetMatrixElement (0, i);
//    //                        dH=matPythag(dF,dG);
//    //                        pWResult->SetMatrixElement (0, i, dH);
//    //                        dH=1.0/dH;
//    //                        dC=dG*dH;
//    //                        dS=(-dF*dH);
//    //                        for (j=0;j<iNumRows;j++) 
//    //                        {
//    //                            dY=a[j][iNM];
//    //                            dZ=a[j][i];
//    //                            this->SetMatrixElement (j, iNM, dY*dC + dZ*dS);
//    //                            this->SetMatrixElement (j, i,   dZ*dC - dY*dS);
//    //                            a[j][iNM]=dY*dC+dZ*dS;
//    //                            a[j][i]=dZ*dC-dY*dS;
//    //                        }
//    //                    }
//    //                }
//    //            }
//    //            dZ= pWResult->GetMatrixElement (0, k);
//    //            if (l == k) 
//    //            {
//    //                if (dZ < 0.0) 
//    //                {
//    //                    pWResult->SetMatrixElement (0, k, -dZ);
//    //                    for (j=0;j<iNumCols;j++) 
//    //                        pVResult->MultCMatrixElement (j, k, -1);
//    //                }
//    //                break;
//    //            }
//    //            if (iIter >= iMaxIterations) 
//    //                return kMatrixNoSVD;
//    //            dX= pWResult->GetMatrixElement (0, l);
//    //            iNM=k-1;
//    //            dY= pWResult->GetMatrixElement (0, iNM);
//    //            dG=pdRV1[iNM];
//    //            dH=pdRV1[k];
//    //            dF=((dY-dZ)*(dY+dZ)+(dG-dH)*(dG+dH))/(2.0*dH*dY);
//    //            dG=matPythag(dF,1.0);
//    //            dF=((dX-dZ)*(dX+dZ)+dH*((dY/(dF+MHLP_SIGN(dG,dF)))-dH))/dX;
//    //            dC=dS=1.0;
//    //            for (j=l;j<=iNM;j++) 
//    //            {
//    //                int jj;
//    //
//    //                i=j+1;
//    //                dG=pdRV1[i];
//    //                dY= pWResult->GetMatrixElement (0, i);
//    //                dH=dS*dG;
//    //                dG=dC*dG;
//    //                dZ=matPythag(dF,dH);
//    //                pdRV1[j]=dZ;
//    //                dC=dF/dZ;
//    //                dS=dH/dZ;
//    //                dF=dX*dC+dG*dS;
//    //                dG=dG*dC-dX*dS;
//    //                dH=dY*dS;
//    //                dY=dY*dC;
//    //                for (jj=0;jj<iNumCols;jj++) 
//    //                {
//    //                    dX= pVResult->GetMatrixElement (jj, j);
//    //                    dZ= pVResult->GetMatrixElement (jj, i);
//    //
//    //                    pVResult->SetMatrixElement (jj, j, dX*dC + dZ*dS);
//    //                    pVResult->SetMatrixElement (jj, i, dZ*dC - dX*dS);
//    //                }
//    //                dZ=matPythag(dF,dH);
//    //                pWResult->SetMatrixElement (0, j, dZ);
//    //                if (dZ) 
//    //                {
//    //                    dZ=1.0/dZ;
//    //                    dC=dF*dZ;
//    //                    dS=dH*dZ;
//    //                }
//    //                dF=(dC*dG)+(dS*dY);
//    //                dX=(dC*dY)-(dS*dG);
//    //                for (jj=0;jj<iNumRows;jj++) 
//    //                {
//    //                    dY=a[jj][j];
//    //                    dZ=a[jj][i];
//    //                    this->SetMatrixElement (jj, j,  dY*dC + dZ*dS);
//    //                    this->SetMatrixElement (jj, i,  dZ*dC - dY*dS);
//    //                    a[jj][j]=dY*dC+dZ*dS;
//    //                    a[jj][i]=dZ*dC-dY*dS;
//    //                }
//    //            }
//    //            pdRV1[l]=0.0;
//    //            pdRV1[k]=dF;
//    //            pWResult->SetMatrixElement (0, k, dX);
//    //
//    //        }
//    //    }
//    //
//    //    delete [] pdRV1;
//
//
//
//
//
//
//
//
//
//
//
//    CMatrix::MatrixError_t rErr;
//    int        i,j,k,l = 0;
//
//    float    *pdRV1; 
//    float    dTmp, 
//        dF,
//        dH, 
//        dS,
//        dNorm       = 0.0,
//        dG          = 0.0,
//        dScale      = 0.0;
//    int        iNumRows   = m_aiMatrixDimensions[kRow],
//        iNumCols   = m_aiMatrixDimensions[kCol];
//
//    if (iNumRows < iNumCols) 
//        return kMatrixNoSVD;
//
//    pdRV1       = new float [iNumCols]; 
//    if (!pdRV1)
//        return kMatrixMemAllocError;
//    rErr    = pWResult->SetDimensions (1, iNumCols);
//    rErr    = pVResult->SetDimensions (iNumCols, iNumCols);
//
//
//    for (i = 0; i < iNumCols; i++) 
//    {
//        l           = i+1;
//        pdRV1[i]    = dScale * dG;
//
//        dG          = 0;
//        dS          = 0;
//        dScale      = 0;
//
//        if (i < iNumRows) 
//        {
//            for (k = i; k < iNumRows; k++) 
//                dScale += fabs (this->GetMatrixElement (k, i));
//
//            if (dScale) 
//            {
//                this->ColMultC (i, i, 1/dScale);
//
//                for (k = i; k < iNumRows; k++) 
//                {
//                    dTmp    = this->GetMatrixElement (k, i);
//                    dS     += dTmp * dTmp;
//                }
//
//                dF      = this->GetMatrixElement (i, i);
//                dG      = -MHLP_SIGN (sqrt (dS), dF);
//                dH      = dF * dG - dS;
//                rErr    = this->SetMatrixElement (i, i, dF-dG);
//
//                if (i != iNumCols) 
//                {
//                    for (j = l; j < iNumCols; j++) 
//                    {
//                        dS      = 0;
//                        for (k=i;k<iNumRows;k++) 
//                            dS     += this->GetMatrixElement (k, i) * this->GetMatrixElement (k, j);
//                        dF      = dS / dH;
//                        for (k=i;k<iNumRows;k++) 
//                            rErr    = this->Add2MatrixElement (k, j, dF * this->GetMatrixElement (k, i));
//                    }
//                }
//                ColMultC (i, i, dScale);
//            }
//        }
//
//        rErr    = pWResult->SetMatrixElement (0, i, dScale*dG);
//
//        dG      = 0;
//        dS      = 0;
//        dScale  = 0;
//
//        if ((i <= iNumRows) && (i != iNumCols)) 
//        {
//            for (k = l; k < iNumCols; k++) 
//                dScale += fabs (this->GetMatrixElement (i, k));
//
//            if (dScale) 
//            {
//                RowMultC (i, l, 1/dScale);
//                for (k=l;k<iNumCols;k++) 
//                {
//                    dTmp    = this->GetMatrixElement (i, k);
//                    dS     += dTmp * dTmp;
//                }
//
//                dF      = this->GetMatrixElement (i, l);
//                dG      = -MHLP_SIGN (sqrt (dS), dF);
//                dH      = dF * dG - dS;
//                rErr    = this->SetMatrixElement (i, l, dF-dG);
//
//                for (k = l; k < iNumCols; k++) 
//                    pdRV1[k]    = this->GetMatrixElement (i, k) / dH;
//
//                if (i != iNumRows) 
//                {
//                    for (j = l; j < iNumRows; j++) 
//                    {
//                        dS      = 0;
//                        for (k = l; k < iNumCols; k++) 
//                            dS     += this->GetMatrixElement (j, k) * this->GetMatrixElement (i, k);
//                        for (k = l; k < iNumCols; k++) 
//                            rErr    = this->Add2MatrixElement (j, k, dS * pdRV1[k]);
//                    }
//                }
//                RowMultC (i, l, dScale);
//            }
//        }
//        dTmp    = fabs (pWResult->GetMatrixElement (0, i)) + fabs (pdRV1[i]);
//        dNorm   = std::max (dNorm, dTmp);
//    }
//
//    for (i = iNumCols-1; i >= 0; i--) 
//    {
//        if (i < iNumCols) 
//        {
//            if (dG) 
//            {
//                for (j = l; j < iNumCols; j++)
//                {
//                    dTmp    = this->GetMatrixElement (i, j) / this->GetMatrixElement (i, l);
//                    rErr    = pVResult->SetMatrixElement (j, i, dTmp / dG);
//                }
//                for (j = l; j < iNumCols; j++)
//                {
//                    dS  = 0;
//                    for (k = l; k < iNumCols; k++) 
//                        dS     += this->GetMatrixElement (i, k) * pVResult->GetMatrixElement (k, j);
//                    for (k = l; k < iNumCols; k++) 
//                        rErr    = pVResult->Add2MatrixElement (k, j, dS * pVResult->GetMatrixElement (k, i));
//                }
//            }
//            for (j = l; j < iNumCols; j++) 
//            {
//                rErr    = pVResult->SetMatrixElement (i, j, 0);
//                rErr    = pVResult->SetMatrixElement (j, i, 0);
//            }
//        }
//        rErr    = pVResult->SetMatrixElement (i, i, 1);
//        dG      = pdRV1[i];
//        l       = i;
//    }
//
//    for (i = iNumCols-1; i >= 0; i--) 
//    {
//        l   = i+1;
//        dG  = pWResult->GetMatrixElement (0, i);
//
//        if (i < iNumCols)
//            RowMultC (i, l, 0);
//
//        if (dG) 
//        {
//            dG  = 1.0 / dG;
//            if (i != iNumCols) 
//            {
//                for (j = l; j < iNumCols; j++) 
//                {
//                    dS  = 0;
//                    for (k = l; k < iNumRows; k++) 
//                        dS     += this->GetMatrixElement (k, i) * this->GetMatrixElement (k, j);
//                    dF  = (dS / this->GetMatrixElement (i, i)) * dG;
//                    for (k = i; k < iNumRows; k++) 
//                        rErr    = this->Add2MatrixElement (k, j, dF * this->GetMatrixElement (k, i));
//                }
//            }
//        } 
//
//        ColMultC (i, i, dG);
//        rErr    = this->Add2MatrixElement (i, i, 1);
//    }
//
//    for (k = iNumCols-1; k >= 0; k--) 
//    {
//        for (int iIter = 0; iIter < iMaxIterations; iIter++) 
//        {
//            float    dC, dX, dY, dZ;
//            int        iNM     = 0,
//                iFlag   = 1;
//            for (l = k; l >= 0; l--) 
//            {
//                iNM = l-1;
//                if (fabs (pdRV1[l]) + dNorm == dNorm) 
//                {
//                    iFlag   = 0;
//                    break;
//                }
//                if (fabs (pWResult->GetMatrixElement (0, iNM)) + dNorm == dNorm) 
//                    break;
//            }
//            if (iFlag) 
//            {
//                dC  = 0.0;
//                dS  = 1.0;
//                for (i = l; i <= k; i++) 
//                {
//                    dF  = dS * pdRV1[i];
//                    if (fabs (dF) + dNorm != dNorm) 
//                    {
//                        dG  = pWResult->GetMatrixElement (0, i);
//                        dH  = matPythag (dF,dG);
//                        pWResult->SetMatrixElement (0, i, dH);
//                        dH  = 1.0 / dH;
//                        dC  = dG * dH;
//                        dS  = (-dF * dH);
//
//                        for (j = 0; j < iNumRows; j++) 
//                        {
//                            dY      = this->GetMatrixElement (j, iNM);
//                            dZ      = this->GetMatrixElement (j, i);
//                            rErr    = this->SetMatrixElement (j, iNM, dY*dC + dZ*dS);
//                            rErr    = this->SetMatrixElement (j, i,   dZ*dC - dY*dS);
//                        }
//                    }
//                }
//            }
//
//            dZ  = pWResult->GetMatrixElement (0, k);
//
//            if (l == k) 
//            {
//                if (dZ < 0.0) 
//                {
//                    rErr    = pWResult->SetMatrixElement (0, k, -dZ);
//                    for (j = 0; j < iNumCols; j++) 
//                        rErr    = pVResult->MultCMatrixElement (j, k, -1);
//                }
//                break;
//            }
//
//            if (iIter >= iMaxIterations) 
//                return kMatrixNoSVD;
//
//            dX  = pWResult->GetMatrixElement (0, l);
//            iNM = k-1;
//            dY  = pWResult->GetMatrixElement (0, iNM);
//            dG  = pdRV1[iNM];
//            dH  = pdRV1[k];
//            dF  = ((dY-dZ) * (dY+dZ) + (dG-dH) * (dG+dH)) / (2.0*dH*dY);
//            dG  = matPythag(dF, 1.0);
//            dF  = ((dX-dZ) * (dX+dZ) + dH*((dY / (dF + MHLP_SIGN(dG,dF)))-dH)) / dX;
//            dC  = 1;
//            dS  = 1;
//
//            for (j = l; j <= iNM; j++) 
//            {
//                int jj;
//
//                i       = j+1;
//                dG      = pdRV1[i];
//                dY      = pWResult->GetMatrixElement (0, i);
//                dH      = dS * dG;
//                dG      = dC * dG;
//                dZ      = matPythag(dF,dH);
//                pdRV1[j]= dZ;
//                dC      = dF / dZ;
//                dS      = dH / dZ;
//                dF      = (dX * dC) + (dG * dS);
//                dG      = (dG * dC) - (dX * dS);
//                dH      = dY * dS;
//                dY      = dY * dC;
//
//                for (jj = 0; jj < iNumCols; jj++) 
//                {
//                    dX      = pVResult->GetMatrixElement (jj, j);
//                    dZ      = pVResult->GetMatrixElement (jj, i);
//
//                    rErr    = pVResult->SetMatrixElement (jj, j, dX*dC + dZ*dS);
//                    rErr    = pVResult->SetMatrixElement (jj, i, dZ*dC - dX*dS);
//                }
//
//                dZ      = matPythag(dF, dH);
//                rErr    = pWResult->SetMatrixElement (0, j, dZ);
//
//                if (dZ) 
//                {
//                    dZ  = 1.0/dZ;
//                    dC  = dF * dZ;
//                    dS  = dH * dZ;
//                }
//                dF  = (dC * dG) + (dS * dY);
//                dX  = (dC * dY) - (dS * dG);
//
//                for (jj = 0; jj < iNumRows; jj++) 
//                {
//                    dY      = this->GetMatrixElement (jj, j);
//                    dZ      = this->GetMatrixElement (jj, i);
//                    rErr    = this->SetMatrixElement (jj, j,  dY*dC + dZ*dS);
//                    rErr    = this->SetMatrixElement (jj, i,  dZ*dC - dY*dS);
//                }
//            }
//            pdRV1[l]    = 0;
//            pdRV1[k]    = dF;
//            rErr        = pWResult->SetMatrixElement (0, k, dX);
//        }
//    }
//
//    delete [] pdRV1;
//
//    return kMatrixNoError;
//}
//
//
//CMatrix::MatrixError_t CMatrix::GetSVD (CMatrix *pUResult, CMatrix *pWResult, CMatrix *pVResult, int iMaxIterations)
//{
//    CMatrix::MatrixError_t rErr;
//
//    rErr    = pUResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pUResult->CalcSVD (pWResult, pVResult, iMaxIterations);
//}
//
////int matrix::QR_householder(matrix &Q, matrix &R){
////    matrix H(*this);
////    // initialize the matrices
////    Q.identity();
////    R=(*this);
////    // find each of the Householder reflection matrices
////    for(int col=0;col<dim_spalten-1;++col){
////        // get the x vector
////        matrix x(dim_spalten - col, 1);
////        x=R.get_submatrix(col, col, dim_spalten-1, col);
////        // make x into u~ (it only involves changing one element, so there's no point in allocating memory for a brand new vector)
////        x.elem[0][0]+=sgn(x.elem[0][0], 1)*x.m2v().norm();
////        // create H_temp (the 2uu^t part of H)
////        matrix H_temp(x.dim_spalten, x.dim_spalten);
////        H_temp=x*x.transpose()*(2.0/((x.transpose()*x).m2s()));
////        // create H (I - 2uu^t)
////        H.identity();
////        householder_corner_subtraction(H, H_temp);
////        // "add" h to the q and r matrices
////        Q*=H;
////        R=H*R;
////        H_temp.deallocate();
////        x.deallocate();
////    }
////    H.deallocate();
////    return 1;
////}
//#ifdef OLDQR
//CMatrix::MatrixError_t   CMatrix::CalcQR (CMatrix *pRResult)
//{
//    CMatrix::MatrixError_t rErr = kMatrixNoError;
//    CMatrix Q, 
//        R, 
//        H, HTmp,
//        X;
//    float *pdTemp    = new float [m_aiMatrixDimensions[kRow]];
//    int    j,
//        iNumOfRows  = m_aiMatrixDimensions[kRow],
//        iNumOfCols  = m_aiMatrixDimensions[kCol];
//
//    //if (iNumOfRows != iNumOfCols)
//    //    return kMatrixNotSquare;
//
//    // initialize the matrices 
//    rErr    = Q.SetDimensions (iNumOfRows, iNumOfRows);
//    rErr    = R.SetDimensions (iNumOfRows, iNumOfCols);
//    rErr    = H.SetDimensions (iNumOfRows, iNumOfRows);
//
//    rErr    = Q.SetMatrixEye ();
//    rErr    = R.Copy (*this);
//
//    //matrix_empty(&x);
//    //matrix_empty(&h_temp);
//
//    // find each of the Householder reflection matrices 
//    for (j = 0; j < iNumOfCols-1; j++)
//    {
//        CMatrix XTmp;
//        // get the x vector 
//        rErr    = X.SetDimensions (iNumOfCols - j, 1);
//        rErr    = R.GetMatrixCol (j, j, pdTemp, iNumOfCols - j);
//        rErr    = X.SetMatrixCol (0, pdTemp, iNumOfCols - j);
//
//        rErr    = XTmp.Copy (X);
//        rErr    = XTmp.CalcTransposition ();
//        rErr    = XTmp.CalcMultiply (&X);
//
//        // make x into u~ (it only involves changing one element, so there's
//        //no point in allocating memory for a brand new vector) 
//        if (X.GetMatrixElement (0, 0) >= 0)
//            X.SetMatrixElement (0, 0, X.GetMatrixElement (0, 0) + sqrt (XTmp.GetMatrixElement (0, 0)));
//        else
//            X.SetMatrixElement (0, 0, X.GetMatrixElement (0, 0) - sqrt (XTmp.GetMatrixElement (0, 0)));
//        //XTmp.MatrixAbs ();
//
//        //X.SetMatrixElement (0, 0, 
//        //    X.GetMatrixElement (0, 0) + XTmp.CalcColMean (0)*XTmp.GetNumRows ()*((X.GetMatrixElement (0, 0) >= 0)? 1 : -1));
//
//        rErr    = XTmp.Copy (X);
//        rErr    = XTmp.CalcTransposition ();
//        rErr    = XTmp.CalcMultiply (&X);
//
//        // create h_temp (the 2uu^t part of h) 
//        //        H_temp=x*x.transpose()*(2.0/((x.transpose()*x).m2s()));
//        rErr    = HTmp.Copy (X);
//        rErr    = X.CalcTransposition ();
//        rErr    = HTmp.CalcMultiply (&X);
//        HTmp.MatrixMultC (2.0 / XTmp.GetMatrixElement (0, 0));
//
//        // create h (I - 2uu^t) 
//        rErr    = H.SetMatrixEye ();
//        rErr    = H.HouseholderCornerSubtraction (HTmp);
//
//        rErr    = Q.CalcMultiply (&H);
//        rErr    = HTmp.Copy (R);
//        rErr    = H.GetMultiply (&R, &HTmp);
//        //rErr    = H.CalcMultiply (&R);
//        //rErr    = R.Copy (H);
//    }
//
//    if (pRResult)
//    {
//        rErr    = pRResult->Copy (R);
//        if (rErr != kMatrixNoError)
//            return rErr;
//    }
//
//    return this->Copy (Q);
//}
//#else
//CMatrix::MatrixError_t   CMatrix::CalcQR (CMatrix *pRResult)
//{
//    //http://matlabdb.mathematik.uni-stuttgart.de/download.jsp?MC_ID=3&MP_ID=164
//    //function [A] = QR_HOUSE(A);
//    //[n,m] = size(A);
//    //for k = 1:min(n-1,m)
//    //    v(k:n,1) = HOUSEHOLDER(A(k:n,k));
//    //    A(k:n,k:m) = HOUSEHOLDER_MULT(A(k:n,k:m),v(k:n,1));
//    //    A(k+1:n,k) = v(k+1:n,1);
//    //end
//
//    CMatrix::MatrixError_t   rErr = kMatrixNoError;
//    float *pdTemp    = new float [m_aiMatrixDimensions[kRow]];
//    int    j,
//        iStop,
//        iNumOfRows  = m_aiMatrixDimensions[kRow],
//        iNumOfCols  = m_aiMatrixDimensions[kCol];
//    CMatrix V,A,R;
//
//    iStop   = std::min(iNumOfRows-1, iNumOfCols);
//    rErr    = R.Copy (*this);
//
//
//    for (j = 0; j < iStop; j++)
//    {
//        rErr    = R.GetMatrixCol (j, j, pdTemp, iNumOfRows - j);
//        A.SetDimensions (iNumOfRows-j,1);
//        rErr    = A.SetMatrixCol (0, pdTemp, iNumOfRows - j);
//        rErr    = V.HouseHolder (A); // restricted size compare to matlab code
//        rErr    = R.GetPartialMatrix (&A, j,j);
//        rErr    = A.HouseHolderMult (V); // restricted size compared to matlab code
//        for (int i = 0; i < iNumOfRows-j; i++)
//        {
//            for (int k = 0; k < iNumOfCols-j; k++)
//                rErr= R.SetMatrixElement (j+i, j+k, A.GetMatrixElement (i,k));
//        }
//        rErr    = V.GetMatrixCol (0, 1, pdTemp, iNumOfRows - j-1);
//        rErr    = R.SetMatrixCol (j, pdTemp, iNumOfRows - j-1, j+1);
//    }
//
//    delete [] pdTemp;
//
//    if (pRResult)
//    {
//        rErr    = pRResult->Copy (R);
//        if (rErr != kMatrixNoError)
//            return rErr;
//    }
//
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::HouseHolder (CMatrix &Input)
//{
//    //function [v] = HOUSEHOLDER(a);
//    //n = length(a);
//    //v = a;
//    //if (a(1) >= 0) beta = a(1) + norm(a);
//    //else beta = a(1) - norm(a);
//    //end
//    //    v(2:n) = 1/beta * v(2:n);
//    //v(1) = 1;
//    CMatrix::MatrixError_t   rErr = kMatrixNoError;
//    float    dMatrixElement = Input.GetMatrixElement (0,0),
//        dNorm;
//    rErr    = this->Copy (Input);
//    rErr    = Input.CalcTransposition ();
//    rErr    = Input.CalcScalarProd (*this, &dNorm);
//
//    if (dNorm >= kSingularityThresh)
//    {
//        dMatrixElement += (dMatrixElement >=0)? sqrt (dNorm) : -sqrt(dNorm);
//        if (abs (dMatrixElement) >= kSingularityThresh)
//            this->MatrixMultC (1.0/dMatrixElement);
//    }
//    this->SetMatrixElement (0, 0, 1.0);
//
//    return kMatrixNoError;
//}
//CMatrix::MatrixError_t CMatrix::HouseHolderMult (CMatrix &VectorInput)
//{
//    //function [A] = HOUSEHOLDER_MULT(A,v);
//    //beta = -2/(v'*v);
//    //w = v'*A;          % w is a line vector
//    //A = A + beta*v*w;
//    CMatrix::MatrixError_t   rErr = kMatrixNoError;
//    float    dBeta;
//    CMatrix     W,Tmp;
//    rErr    = W.Copy (VectorInput);
//    rErr    = W.CalcTransposition ();
//    rErr    = W.CalcScalarProd (VectorInput, &dBeta);
//    dBeta = -2.0/dBeta;
//    rErr    = W.CalcMultiply (this);
//    rErr    = VectorInput.GetMultiply (&Tmp, &W);
//    Tmp.MatrixMultC (dBeta);
//    rErr    = this->CalcAdd (&Tmp);
//
//    return kMatrixNoError;
//}
//CMatrix::MatrixError_t CMatrix::GetPartialMatrix (CMatrix *pResult, int iStartRow, int iStartCol)
//{
//    int i,i2;
//    if ((iStartRow >= m_aiMatrixDimensions[kRow]) ||(iStartCol >= m_aiMatrixDimensions[kCol]))
//        return kMatrixDimensionMisfit;
//
//    pResult->SetDimensions (m_aiMatrixDimensions[kRow]-iStartRow, m_aiMatrixDimensions[kCol]-iStartCol);
//
//    for (i = iStartRow,i2 = 0; i < m_aiMatrixDimensions[kRow]; i++,i2++)
//        pResult->SetMatrixRow (i2, &m_ppfMatrix[i][iStartCol], m_aiMatrixDimensions[kCol]-iStartCol);
//
//    return kMatrixNoError;
//}
//#endif
//CMatrix::MatrixError_t CMatrix::CalcScalarProd (CMatrix &V, float *pdResult)
//{
//    int i;
//    if (!pdResult)
//        return kMatrixIllegalPointer;
//    if (m_aiMatrixDimensions[kRow] != 1 || V.GetNumColumns () != 1 || m_aiMatrixDimensions[kCol] != V.GetNumRows ())
//        return kMatrixDimensionMisfit;
//
//    *pdResult = 0;
//
//    for (i = 0; i < m_aiMatrixDimensions[kCol]; i++)
//        *pdResult  += m_ppfMatrix[0][i] * V.GetMatrixElement(i,0);
//
//    return kMatrixNoError;
//}
//CMatrix::MatrixError_t CMatrix::GetQR (CMatrix *pQResult, CMatrix *pRResult)
//{
//    CMatrix::MatrixError_t rErr;
//
//    rErr    = pQResult->Copy(*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pQResult->CalcQR (pRResult);
//}
//
//CMatrix::MatrixError_t   CMatrix::CalcInversion ()
//{
//    CMatrix::MatrixError_t rErr;
//    CMatrix TmpMatrix,
//        EyeMatrix;
//    int    i,j,
//        iNumRows    = m_aiMatrixDimensions[kRow],
//        iNumCols    = m_aiMatrixDimensions[kCol];
//    float dDet       = 1;
//
//    if (iNumRows!=iNumCols)
//        return kMatrixNotSquare;
//
//    rErr    = TmpMatrix.Copy (*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    rErr    = EyeMatrix.SetDimensions (iNumRows, iNumCols);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    EyeMatrix.SetMatrixEye ();
//
//    if (m_ppfMatrix[0][0] == 0)
//    {
//        i   =1;
//        while (i < iNumRows)
//        {
//            if (m_ppfMatrix[i][0]!=0)
//            {
//                this->SwapRowCol(0,i);
//                EyeMatrix.SwapRowCol(0,i);
//                dDet *= -1;
//                break;
//            }
//            i++;
//        }			
//    }
//
//    dDet   *= m_ppfMatrix[0][0];
//
//    EyeMatrix.RowDivC (0, m_ppfMatrix[0][0]);
//    this->RowDivC (0, m_ppfMatrix[0][0]);
//
//    for (i = 1; i < iNumRows; i++)
//    {
//        j   = 0;
//
//        if (dDet < kSingularityThresh*kSingularityThresh)
//            dDet = 0;
//
//        while (j < i)
//        {
//            EyeMatrix.RowSub (i, j, m_ppfMatrix[i][j]);
//            this->RowSub (i, j, m_ppfMatrix[i][j]);
//            j++;
//        }
//
//        if (m_ppfMatrix[i][i] != 0)
//        {
//            dDet   *= m_ppfMatrix[i][i];
//            EyeMatrix.RowDivC (i, m_ppfMatrix[i][i]); 
//            this->RowDivC (i, m_ppfMatrix[i][i]); 
//
//        }
//
//        if (m_ppfMatrix[i][i] == 0)
//        {
//            for (int j1 = i+1; j1 < iNumCols; j1++)
//            {
//                if (m_ppfMatrix[i][j1] != 0)			// Column pivotting not supported
//                {
//                    for (int i1=0; i1 < iNumRows; i1++)
//                    {
//                        for (j = 0; j < iNumCols; j++)
//                            m_ppfMatrix[i1][j]  = TmpMatrix.GetMatrixElement (i1, j);
//                    }
//                    rErr = TmpMatrix.Inverse ();
//                    if (rErr == kMatrixNoError)
//                        this->Copy (TmpMatrix);
//
//                    return rErr;
//                }
//            }
//        }
//    }
//
//    for (i = iNumRows-1; i > 0; i--)
//    {
//        for (j = i-1; j >= 0; j--)
//        {
//            EyeMatrix.RowSub (j, i, m_ppfMatrix[j][i]);
//            this->RowSub (j, i, m_ppfMatrix[j][i]);
//        }
//    }						
//
//    return this->Copy (EyeMatrix);
//}
//
//CMatrix::MatrixError_t CMatrix::GetInverted (CMatrix *pResult)
//{
//    CMatrix::MatrixError_t rErr;
//
//    rErr    = pResult->Copy (*this);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return pResult->CalcInversion ();
//}
//
//float CMatrix::CalcTrace ()
//{
//    int        iStop   = std::min (m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);
//    float    dResult = 0.0;
//
//    for (int i = 0; i < iStop; i++)
//        dResult    += m_ppfMatrix[i][i];
//
//    return dResult;
//}
//
//float CMatrix::CalcDeterminant ()
//{
//    CMatrix::MatrixError_t rErr;
//    CMatrix TmpMatrix;
//    int    i,j,
//        iNumRows    = m_aiMatrixDimensions[kRow],
//        iNumCols    = m_aiMatrixDimensions[kCol];
//    float dDet       = 1;
//
//    rErr    = TmpMatrix.Copy (*this);
//
//    if (iNumRows == 2)
//        return( (m_ppfMatrix[0][0] * m_ppfMatrix[1][1]) - (m_ppfMatrix[0][1] * m_ppfMatrix[1][0])); 
//
//    if (m_ppfMatrix[0][0] == 0)
//    {
//        i   = 1;
//        while (i<iNumRows)
//        {
//            if (m_ppfMatrix[i][0] != 0)
//            {
//                this->SwapRowCol(0,i);
//                dDet *= -1;
//                break;
//            }
//            i++;
//        }
//    }
//
//    if (m_ppfMatrix[0][0] == 0)
//    {
//        rErr    = this->Copy (TmpMatrix);
//
//        return 0;							//If all the elements in a row or column of matrix are 0, determient is equal to 0
//    }
//
//    dDet   *= m_ppfMatrix[0][0];
//
//    this->RowDivC (0, m_ppfMatrix[0][0]);
//
//    for (i = 1; i < iNumRows; i++)
//    {
//        j=0;
//
//        if (dDet < kSingularityThresh*kSingularityThresh)
//            dDet = 0;
//
//        while (j < i)							//preparing an upper triangular matrix
//        {
//            this->RowAdd (i, j, -m_ppfMatrix[i][j]);
//            j++;
//        }
//
//        if (m_ppfMatrix[i][i] != 0)
//        {
//            dDet   *= m_ppfMatrix[i][i];					//Dividing the entire row with non zero diagonal element. Multiplying det with that factor.	
//            this->RowDivC (i, m_ppfMatrix[i][i]); 
//        }
//
//        if (m_ppfMatrix[i][i] == 0)						// Chcek if the diagonal elements are zeros
//        {
//            for (j = i+1; j < iNumCols; j++)
//            {
//                if (m_ppfMatrix[i][j]!= 0)
//                {
//                    this->ColAdd (i, j, 1);//Adding of columns do not change the determinant
//
//                    dDet *= m_ppfMatrix[i][i];
//                    this->RowDivC (i, m_ppfMatrix[i][i]);
//                    break;
//                }
//            }
//        }
//
//        if (m_ppfMatrix[i][i]== 0)						//if diagonal element is still zero, Determinant is zero.
//        {
//            rErr    = this->Copy (TmpMatrix);
//
//            return 0;							//If all the elements in a row or column of matrix are 0, determient is equal to 0
//        }
//    }
//
//    rErr    = this->Copy (TmpMatrix);
//
//    return dDet;
//}
//
//float CMatrix::CalcRowMean (int iRow)
//{
//    float dResult = 0;
//
//    if (m_aiMatrixDimensions[kCol] <= 0 || m_aiMatrixDimensions[kRow] <= 0)
//        return 0;
//
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//        dResult += m_ppfMatrix[iRow][j];
//
//    return dResult / m_aiMatrixDimensions[kCol];
//}
//
//float CMatrix::CalcColMean (int iCol)
//{
//    float dResult = 0;
//
//    if (m_aiMatrixDimensions[kCol] <= 0 || m_aiMatrixDimensions[kRow] <= 0)
//        return 0;
//
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        dResult += m_ppfMatrix[i][iCol];
//
//    return dResult / m_aiMatrixDimensions[kRow];
//}
//
//bool   CMatrix::IsMatrixSingular ()
//{
//    return (this->CalcDeterminant () < kSingularityThresh);
//}
//
//CMatrix::MatrixError_t CMatrix::SetMatrixZero ()
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//            m_ppfMatrix[i][j]   = 0;
//
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t  CMatrix::SetMatrixEye ()
//{
//    int    iStop = std::min (m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);
//
//    this->SetMatrixZero ();
//    for (int i = 0; i < iStop; i++)
//        m_ppfMatrix[i][i]   = 1;
//
//    return kMatrixNoError;
//}
//
//void   CMatrix::MatrixMultC (float fValue)
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//            m_ppfMatrix[i][j]   *= fValue;
//
//    return;
//}
//void   CMatrix::MatrixAbs ()
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//            m_ppfMatrix[i][j]    = abs (m_ppfMatrix[i][j]);
//
//    return;
//}
//void   CMatrix::MatrixLog ()
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//        {
//            if (m_ppfMatrix[i][j] < kSingularityThresh)
//                m_ppfMatrix[i][j] = kSingularityThresh;
//            m_ppfMatrix[i][j]    = log (m_ppfMatrix[i][j]);
//        }
//
//        return;
//}
//
//void   CMatrix::RowDivC (int iRow, float fValue)
//{
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//        m_ppfMatrix[iRow][j]   /= fValue;
//
//    return;
//}
//
//void   CMatrix::RowMultC (int iRow, float fValue)
//{
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//        m_ppfMatrix[iRow][j]   *= fValue;
//
//    return;
//}
//
//void   CMatrix::RowMultC (int iRow, int iColBegin, float fValue)
//{
//    for (int j = iColBegin; j < m_aiMatrixDimensions[kCol]; j++)
//        m_ppfMatrix[iRow][j]   *= fValue;
//
//    return;
//}
//
//void   CMatrix::RowAdd (int iRow1, int iRow2, float fValue)
//{
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//        m_ppfMatrix[iRow1][j]   += (fValue *m_ppfMatrix[iRow2][j]);
//
//    return;
//}
//
//void   CMatrix::RowSub (int iRow1, int iRow2, float fValue)
//{
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//        m_ppfMatrix[iRow1][j]   -= (fValue *m_ppfMatrix[iRow2][j]);
//
//    return;
//}
//
//void   CMatrix::ColDivC (int iCol, float fValue)
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kCol]; i++)
//        m_ppfMatrix[i][iCol]   /= fValue;
//
//    return;
//}
//
//void   CMatrix::ColMultC (int iCol, float fValue)
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        m_ppfMatrix[i][iCol]   *= fValue;
//
//    return;
//}
//
//void   CMatrix::ColMultC (int iCol, int iRowBegin, float fValue)
//{
//    for (int i = iRowBegin; i < m_aiMatrixDimensions[kRow]; i++)
//        m_ppfMatrix[i][iCol]   *= fValue;
//
//    return;
//}
//
//void   CMatrix::ColAdd (int iCol1, int iCol2, float fValue)
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kCol]; i++)
//        m_ppfMatrix[i][iCol1]   += (fValue * m_ppfMatrix[i][iCol2]);
//
//    return;
//}
//
//void   CMatrix::ColAddC (int iCol, float fValue)
//{
//    for (int i = 0; i < m_aiMatrixDimensions[kCol]; i++)
//        m_ppfMatrix[i][iCol]   += fValue;
//
//    return;
//}
//
//void   CMatrix::SwapRowCol (int iRow, int iCol)
//{
//    for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//    {
//        float dTmp       = m_ppfMatrix[iRow][j];
//        m_ppfMatrix[iRow][j]= m_ppfMatrix[iCol][j];
//        m_ppfMatrix[iCol][j]= dTmp;
//    }
//}
//
//CMatrix::MatrixError_t CMatrix::Inverse ()
//{
//    int    i,j,
//        iNumRows    = m_aiMatrixDimensions[kRow],
//        iNumCols    = m_aiMatrixDimensions[kCol];
//    CMatrix Inverted;
//    CMatrix::MatrixError_t rErr;
//    float dDet = this->CalcDeterminant ();
//
//    if (this->IsMatrixSingular ())
//        return kMatrixIsSingular;
//
//    rErr    = Inverted.SetDimensions (iNumRows, iNumCols);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    rErr    = Inverted.SetMatrixZero ();
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    for (i = 0; i <iNumRows; i++)
//    {
//        for (j = 0; j < iNumCols; j++)
//        {
//            CMatrix SMatrix;
//            float dDetS,
//                dFac = (((i+j+2)%2) == 0)? 1 : -1;
//            rErr    = this->GetSubMatrix (&SMatrix, i, j);
//            if (rErr != kMatrixNoError)
//                return rErr;
//            dDetS   = SMatrix.CalcDeterminant ();
//            rErr    = Inverted.SetMatrixElement (i, j, dFac * dDetS / dDet);
//            if (rErr != kMatrixNoError)
//                return rErr;
//        }				
//    }
//    rErr    = Inverted.CalcTransposition ();
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    rErr    = this->Copy (Inverted);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::GetSubMatrix (CMatrix *pResult, int iRow, int iCol)
//{
//    CMatrix::MatrixError_t rErr;
//    int    iNewSize    = m_aiMatrixDimensions[kRow] - 1;
//    int    i           = 0;
//
//    rErr    = pResult->SetDimensions (iNewSize, iNewSize);
//    if (rErr != kMatrixNoError)
//        return rErr;
//
//    for (int p = 0; p < iNewSize; p++)
//    {				
//        if (p != iRow)
//        {
//            int j  = 0;
//            for(int q = 0; q < iNewSize; q++)
//            {
//                if (q != iCol)
//                {
//                    pResult->SetMatrixElement (i, j, m_ppfMatrix[p][q]);
//                    j  += 1;
//                }
//            }
//            i  += 1;
//        }
//    }
//
//    return kMatrixNoError;
//}
//
////void matrix::householder_corner_subtraction(matrix &left_large, const matrix &right_small){
////    /* calculate start offsets for the large matrix */
////    int row_offset=left_large.dim_zeilen-right_small.dim_zeilen;
////    int col_offset=left_large.dim_spalten-right_small.dim_spalten;
////    /* subtract the elements individually */
////    for(int row=0;row<right_small.dim_zeilen;++row)
////        for(int col=0;col<right_small.dim_spalten;++col)
////            left_large.elem[row_offset+row][col_offset+col]-= right_small.elem[row][col];
////}
//
//CMatrix::MatrixError_t  CMatrix::HouseholderCornerSubtraction (CMatrix &SmallerMatrix)
//{
//    int    iNumSmallRows   = SmallerMatrix.GetNumRows (),
//        iNumSmallCols   = SmallerMatrix.GetNumColumns ();
//    int    i, j,
//        iRowOffset      = m_aiMatrixDimensions[kRow] - iNumSmallRows,
//        iColOffset      = m_aiMatrixDimensions[kCol] - iNumSmallCols;
//
//    if ((m_aiMatrixDimensions[kRow] < iNumSmallRows) || (m_aiMatrixDimensions[kCol] < iNumSmallCols))
//        return kMatrixDimensionMisfit;
//
//
//    // subtract the elements individually 
//    for (i = 0; i < iNumSmallRows; i++)
//    {
//        for (j = 0; j < iNumSmallCols; j++)
//            this->SetMatrixElement(iRowOffset + i, iColOffset + j, this->GetMatrixElement (iRowOffset + i, iColOffset + j) - SmallerMatrix.GetMatrixElement (i, j));
//    }
//
//    return kMatrixNoError;
//}
//
//CMatrix::MatrixError_t CMatrix::WriteMatrix2File (std::string strFileName)
//{
//    std::ofstream aFFileHandle;
//
//    if (strFileName.empty ())
//        return kMatrixNoError;
//
//    aFFileHandle.open (strFileName.c_str (), std::ios::out);
//
//    if (aFFileHandle.is_open ())
//    {
//        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//        {
//            for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
//                aFFileHandle << std::scientific << std::setprecision(10) << this->GetMatrixElement (i, j) << "\t";
//            aFFileHandle << std::endl;
//        }
//    }
//    if (aFFileHandle.is_open ())
//        aFFileHandle.close ();
//
//    return kMatrixNoError;
//}
//
//
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//// operators
//
//
////CMatrix& CMatrix::operator=(CMatrix &Matrix)
////{
////    int    i,
////            iNewNumRows = Matrix.GetNumRows (),
////            iNewNumCols = Matrix.GetNumColumns ();
////
////    CMatrix::MatrixError_t   rErr;
////
////    if (&Matrix == this)
////        return *this;
////
////    //this->Delete ();
////
////    rErr = this->SetDimensions (iNewNumRows, iNewNumCols);
////
////    ZASSERT (rErr != kMatrixNoError);
////
////    for (i = 0; i < iNewNumRows; i++)       
////        Matrix.GetMatrixRow (i, m_ppfMatrix[i], iNewNumCols);
////
////    return *this ;
////}
////
////CMatrix	CMatrix::operator+(CMatrix &Matrix)
////{
////    // construct the object we are going to return
////    CMatrix		MatrixResult;
////
////    // now add in the other matrix
////    int    iNumRows = Matrix.GetNumRows (),
////            iNumCols = Matrix.GetNumColumns ();
////    // first check for a valid addition operation
////    ZASSERT ((m_aiMatrixDimensions[kCol] != iNumCols) || (m_aiMatrixDimensions[kRow] != iNumRows));
////
////    MatrixResult.SetDimensions (iNumRows, iNumCols);
////
////    for (int i = 0 ; i < iNumRows ; ++i)
////    {
////        for (int j = 0 ; j < iNumCols ; ++j)
////            MatrixResult.SetMatrixElement(i, j, this->GetMatrixElement(i, j) + Matrix.GetMatrixElement(i, j)) ;
////    }
////    return MatrixResult ;
////}
////
////CMatrix	CMatrix::operator-(CMatrix &Matrix)
////{
////    // construct the object we are going to return
////    CMatrix		MatrixResult;
////
////    // now add in the other matrix
////    int    iNumRows = Matrix.GetNumRows (),
////            iNumCols = Matrix.GetNumColumns ();
////    // first check for a valid addition operation
////    ZASSERT ((m_aiMatrixDimensions[kCol] != iNumCols) || (m_aiMatrixDimensions[kRow] != iNumRows));
////
////    MatrixResult.SetDimensions (iNumRows, iNumCols);
////
////    for (int i = 0 ; i < iNumRows ; ++i)
////    {
////        for (int j = 0 ; j < iNumCols ; ++j)
////            MatrixResult.SetMatrixElement(i, j, this->GetMatrixElement(i, j) - Matrix.GetMatrixElement(i, j)) ;
////    }
////    return MatrixResult ;
////}
////
////CMatrix	CMatrix::operator*(CMatrix &Matrix)
////{
////    float    fValue;
////    int        iNumNewColumns  = Matrix.GetNumColumns (),
////                iNumOldColumns  = Matrix.GetNumRows ();
////    // construct the object we are going to return
////    CMatrix		MatrixResult;
////
////    ZASSERT (m_aiMatrixDimensions[kCol] != iNumOldColumns);
////
////    MatrixResult.SetDimensions (m_aiMatrixDimensions[kRow], iNumNewColumns);
////
////    // e.g.
////    // [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
////    // [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
////    //             [K][L]
////    //
////    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
////    {
////        for (int j = 0; j < iNumNewColumns; j++)
////        {
////            fValue  = 0.0;
////            for (int k = 0; k < iNumOldColumns ; k++)
////                fValue += GetMatrixElement(i, k) * Matrix.GetMatrixElement(k, j) ;
////
////            MatrixResult.SetMatrixElement(i, j, fValue) ;
////        }
////    }
////    return MatrixResult;
////}
////
////void   CMatrix::operator*=(CMatrix &Matrix)
////{
////    *this   = *this * Matrix;
////}
////
////void	CMatrix::operator+=(CMatrix &Matrix)
////{
////    // now add in the other matrix
////    int    iNumRows = Matrix.GetNumRows (),
////            iNumCols = Matrix.GetNumColumns ();
////    // first check for a valid addition operation
////    ZASSERT ((m_aiMatrixDimensions[kCol] != iNumCols) || (m_aiMatrixDimensions[kRow] != iNumRows));
////
////    for (int i = 0 ; i < iNumRows ; ++i)
////    {
////        for (int j = 0 ; j < iNumCols ; ++j)
////            m_ppfMatrix[i][j]  += Matrix.GetMatrixElement(i, j);
////    }
////    return;
////}
////
////void	CMatrix::operator-=(CMatrix &Matrix)
////{
////    // now add in the other matrix
////    int    iNumRows = Matrix.GetNumRows (),
////            iNumCols = Matrix.GetNumColumns ();
////    // first check for a valid addition operation
////    ZASSERT ((m_aiMatrixDimensions[kCol] != iNumCols) || (m_aiMatrixDimensions[kRow] != iNumRows));
////
////    for (int i = 0 ; i < iNumRows ; ++i)
////    {
////        for (int j = 0 ; j < iNumCols ; ++j)
////            m_ppfMatrix[i][j]  -= Matrix.GetMatrixElement(i, j);
////    }
////    return;
////}
////
////void	CMatrix::operator*=(float fValue)
////{
////    for (int i = 0 ; i < m_aiMatrixDimensions[kRow] ; ++i)
////    {
////        for (int j = 0 ; j < m_aiMatrixDimensions[kCol] ; ++j)
////            m_ppfMatrix[i][j]  *= fValue;
////    }
////    return;
////}
////
////void	CMatrix::operator+=(float fValue)
////{
////    for (int i = 0 ; i < m_aiMatrixDimensions[kRow] ; ++i)
////    {
////        for (int j = 0 ; j < m_aiMatrixDimensions[kCol] ; ++j)
////            m_ppfMatrix[i][j]  += fValue;
////    }
////    return;
////}
//
//bool   CMatrix::operator==(CMatrix &Matrix)
//{
//    float    *pdTempValues = 0;
//    if (&Matrix == this)
//        return true;
//    if (m_aiMatrixDimensions[kRow] != Matrix.GetNumRows ())
//        return false;
//    if (m_aiMatrixDimensions[kCol] != Matrix.GetNumColumns ())
//        return false;
//
//    pdTempValues    = new float [m_aiMatrixDimensions[kCol]];
//
//    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
//    {
//        CMatrix::MatrixError_t   rErr;
//        rErr    = Matrix.GetMatrixRow (i, pdTempValues, m_aiMatrixDimensions[kCol]);
//        if (rErr)
//        {
//            delete [] pdTempValues;
//            return false;
//        }
//        if (memcmp (m_ppfMatrix[i], pdTempValues, m_aiMatrixDimensions[kCol] * sizeof(float)) != 0)
//        {
//            delete [] pdTempValues;
//            return false;
//        }
//    }
//    delete [] pdTempValues;
//    return true;
//
//}
//
//

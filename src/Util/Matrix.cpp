#include <string>
#include <cassert>

#include "Util.h"
#include "SignalGen.h"
#include "Matrix.h"



CMatrix::CMatrix () :
    m_ppfMatrix(0)
{
    m_aiMatrixDimensions[kRow]  = 0;
    m_aiMatrixDimensions[kCol]  = 0;
}

CMatrix::CMatrix( const CMatrix &other ) :
    m_ppfMatrix(0)
{
    assign(other);
}

CMatrix::CMatrix( int iNumRows, int iNumColumns ) :
    m_ppfMatrix(0)
{
    m_aiMatrixDimensions[kRow]  = 0;
    m_aiMatrixDimensions[kCol]  = 0;

    init(iNumRows, iNumColumns);
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


CMatrix::MatrixError_t   CMatrix::init (int iNumRows, int iNumCols)
{
    int i;

    assert (iNumRows >= 0);
    assert (iNumCols >= 0);

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
    if(!isIndexValid(iRow, 0) || iNumOfValues != m_aiMatrixDimensions[kCol] || !pfValues)
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);
    assert (m_ppfMatrix[iRow] != 0);

    CUtil::copyBuff(m_ppfMatrix[iRow], pfValues, iNumOfValues);

    return kMatrixNoError;
}

CMatrix::MatrixError_t   CMatrix::getRow (int iRow, float *pfValues, int iNumOfValues) const
{
    if(!isIndexValid(iRow, 0) || iNumOfValues != m_aiMatrixDimensions[kCol] || !pfValues)
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);
    assert (m_ppfMatrix[iRow] != 0);

    CUtil::copyBuff(pfValues, m_ppfMatrix[iRow], iNumOfValues);

    return kMatrixNoError;
}

const float* CMatrix::getRow( int iRow ) const
{
    if(!isIndexValid(iRow, 0))
        return 0;

    return m_ppfMatrix[iRow];
}

CMatrix::MatrixError_t   CMatrix::setCol (int iCol, const float *pfValues, int iNumOfValues)
{
    if(!isIndexValid(0, iCol) || iNumOfValues != m_aiMatrixDimensions[kRow] || !pfValues)
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);

    for (int i = 0; i < iNumOfValues; i++)
        m_ppfMatrix[i][iCol]    = pfValues[i];

    return kMatrixNoError;
}

CMatrix::MatrixError_t CMatrix::setZero()
{
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::setZero(m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}


CMatrix::MatrixError_t   CMatrix::getCol (int iCol, float *pfValues, int iNumOfValues) const
{
    if(!isIndexValid(0, iCol) || iNumOfValues != m_aiMatrixDimensions[kRow] || !pfValues)
        return kMatrixIllegalFunctionParam;

    assert (m_ppfMatrix != 0);

    for (int i = 0; i < iNumOfValues; i++)
        pfValues[i] = m_ppfMatrix[i][iCol];

    return kMatrixNoError;
}

CMatrix& CMatrix::operator=(const CMatrix &other)
{ 
    if(this == &other)
        return *this;

    init(other.getNumRows(), other.getNumCols());
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        other.getRow(i, m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);
    return *this;
}

CMatrix CMatrix::operator+(const CMatrix &other) const
{ 
    assert(m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols());

    CMatrix TmpMatrix;

    if (m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols())
    { 
        TmpMatrix.init(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
            {
                TmpMatrix.setElement(i, j, m_ppfMatrix[i][j] + other.getElement(i,j));
            }
        }
    }
    return TmpMatrix;
}

CMatrix CMatrix::operator-(const CMatrix &other) const
{ 
    assert(m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols());
    CMatrix TmpMatrix;

    if (m_aiMatrixDimensions[kRow] == other.getNumRows() && m_aiMatrixDimensions[kCol] == other.getNumCols())
    { 
        TmpMatrix.init(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
            {
                TmpMatrix.setElement(i, j, m_ppfMatrix[i][j] - other.getElement(i,j));
            }
        }
    }
    return TmpMatrix;
}

CMatrix CMatrix::operator*(const CMatrix &other) const
{ 
    assert(m_aiMatrixDimensions[kCol] == other.getNumRows());
    CMatrix TmpMatrix;

    if (m_aiMatrixDimensions[kCol] == other.getNumRows())
    { 
        TmpMatrix.init(m_aiMatrixDimensions[kRow], other.getNumCols());

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            for (int j = 0; j < other.getNumCols (); j++)
            {
                float fResult = 0;
                for (int k = 0; k < m_aiMatrixDimensions[kCol]; k++)
                    fResult     += m_ppfMatrix[i][k] * other.getElement(k,j);
                TmpMatrix.setElement(i, j, fResult);
            }
        }
    }
    
    return TmpMatrix;
}

bool   CMatrix::operator==(const CMatrix &other) const
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

CMatrix CMatrix::operator+(const float fValue) const
{ 
    CMatrix TmpMatrix;

    TmpMatrix.init(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
    {
        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
        {
            TmpMatrix.setElement(i, j, m_ppfMatrix[i][j] + fValue);
        }
    }
    return TmpMatrix;
}

CMatrix CMatrix::operator-(const float fValue) const
{ 
    CMatrix TmpMatrix;

    TmpMatrix.init(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
    {
        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
        {
            TmpMatrix.setElement(i, j, m_ppfMatrix[i][j] - fValue);
        }
    }
    return TmpMatrix;
}

CMatrix CMatrix::operator*(const float fValue) const
{ 
    CMatrix TmpMatrix;
    TmpMatrix.init(m_aiMatrixDimensions[kRow], m_aiMatrixDimensions[kCol]);

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
    {
        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
        {
            TmpMatrix.setElement(i, j, m_ppfMatrix[i][j] * fValue);
        }
    }
    return TmpMatrix;
}

CMatrix::MatrixError_t CMatrix::transpose_I()
{
    this->assign(this->transpose());

    return kMatrixNoError;
}

CMatrix CMatrix::transpose()
{
    int    iNumOldRows = m_aiMatrixDimensions[kRow],
        iNumOldCols = m_aiMatrixDimensions[kCol],
        iNumNewRows = m_aiMatrixDimensions[kCol],
        iNumNewCols = m_aiMatrixDimensions[kRow];
    CMatrix NewMatrix;
    NewMatrix.init(iNumNewRows, iNumNewCols);

    for (int i = 0; i < iNumOldRows; i++)
    {
        for (int j = 0; j < iNumOldCols; j++)
        {
            MatrixError_t rErr = NewMatrix.setElement(j, i, this->getElement(i, j));

            assert (rErr == kMatrixNoError);
            //if (rErr != kMatrixNoError)
            //    return rErr;
        }
    }

    return NewMatrix;
}

void CMatrix::assign( const CMatrix &other )
{
    init(other.getNumRows(), other.getNumCols());
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        other.getRow(i, m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);
}

CMatrix::MatrixError_t CMatrix::mulByElement_I( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::mulBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}


CMatrix::MatrixError_t CMatrix::divByElement_I( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::divBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

CMatrix::MatrixError_t CMatrix::addByElement_I( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::addBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

CMatrix::MatrixError_t CMatrix::subByElement_I( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::subBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

CMatrix CMatrix::mulByElement( const CMatrix &other ) const
{
    CMatrix Result = CMatrix(*this);

    Result.mulByElement_I(other);
    return Result;
}

CMatrix CMatrix::divByElement( const CMatrix &other ) const
{
    CMatrix Result = CMatrix(*this);

    Result.divByElement_I(other);
    return Result;
}

CMatrix CMatrix::addByElement( const CMatrix &other ) const
{
    CMatrix Result = CMatrix(*this);

    Result.addByElement_I(other);
    return Result;
}

CMatrix CMatrix::subByElement( const CMatrix &other ) const
{
    CMatrix Result = CMatrix(*this);

    Result.subByElement_I(other);
    return Result;
}

float CMatrix::getNorm( int p /*= 1*/ ) const
{
    if (p<=0)
        return -1;

    float fResult = 0;

    if (p==1)
    {
        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
        {
            float fTmp = getVectorNorm(-1,j,p);
            //for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
            //    fTmp += abs(m_ppfMatrix[i][j]);

            if (fTmp > fResult)
                fResult = fTmp;
        }
    }
    else if (p == 2)
    {
        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            //float fTmp  = getVectorNorm(i,-1,p);
            //fResult    += fTmp * fTmp;
            fResult += CUtil::mulBuffScalar(m_ppfMatrix[i], m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);
        }

        fResult = sqrt(fResult);
    }
    else
    {
        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
            for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
                fResult += pow(m_ppfMatrix[i][j], 1.F*p);

        fResult = pow(fResult, 1.F/p);
    }

    return fResult;
}

CMatrix::MatrixError_t CMatrix::normalize_I( ActionAppliedTo_t eActionArea /*= kAll*/, int p /*= 1*/ )
{
    float fNorm;
    
    if (eActionArea == kAll)
    {
        fNorm   = getNorm(p);

        if (fNorm <= 0)
            return kMatrixNoError;

        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
            CUtil::mulBuffC(m_ppfMatrix[i], 1.F/fNorm, m_aiMatrixDimensions[kCol]);
    }
    else if (eActionArea == kRow)
    {
        for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        {
            fNorm   = getVectorNorm(i,-1,p);
            if (fNorm <= 0)
                return kMatrixNoError;

            CUtil::mulBuffC(m_ppfMatrix[i], 1.F/fNorm, m_aiMatrixDimensions[kCol]);
        }
    }
    else if (eActionArea == kCol)
    {
        for (int j = 0; j < m_aiMatrixDimensions[kCol]; j++)
        {
            fNorm   = getVectorNorm(-1,j,p);
            if (fNorm <= 0)
                return kMatrixNoError;

            for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
                m_ppfMatrix[i][j]  /= fNorm;
        }
    }

    return kMatrixNoError;
}

CMatrix CMatrix::normalize( ActionAppliedTo_t eActionArea /*= kAll*/, int p /*= 1*/ )
{
    CMatrix Result = *this;
    Result.normalize_I(eActionArea, p);

    return Result;
}

CMatrix::MatrixError_t CMatrix::setRand()
{
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CSignalGen::generateNoise(m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

float CMatrix::getSum( bool bAbs /*= false*/ ) const
{
    float fResult = 0;
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        fResult += CUtil::sum(m_ppfMatrix[i], m_aiMatrixDimensions[kCol], bAbs);

    return fResult;
}

CMatrix CMatrix::mulByOnes( int iNumRows, int iNumCols )
{
    CMatrix Result;
    int iThisRows = getNumRows();

    if (getNumCols() != iNumRows)
        return Result;

    Result.init(iThisRows, iNumCols);

    for (int i = 0; i < iThisRows; i++)
    {
        float fSum = CUtil::sum(m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);
        CUtil::setValue(Result.m_ppfMatrix[i], fSum, iNumCols);
    }
    return Result;
}
CMatrix::MatrixError_t   CMatrix::mulByOnes_I(int iNumRows, int iNumCols)
{
    CMatrix Result = this->mulByOnes(iNumRows, iNumCols);
    this->assign(Result);

    return kMatrixNoError;
}

float CMatrix::getVectorNorm( int iRow /*= -1*/, int iCol /*= -1*/, int p /*= 1*/ ) const
{
    float fResult = 0;
    assert ((isIndexValid(0, iCol) || isIndexValid(iRow, 0)) || (p > 0));

    if (p == 1)
    {
        if (iRow >= 0)
        {
            fResult = CUtil::sum(m_ppfMatrix[iRow], m_aiMatrixDimensions[kCol], true);
        }
        else if (iCol >= 0)
        {
            for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
                fResult += abs(m_ppfMatrix[i][iCol]);
        }
    }
    else if (p == 2)
    {
        if (iRow >= 0)
        {
            fResult = CUtil::mulBuffScalar(m_ppfMatrix[iRow], m_ppfMatrix[iRow], m_aiMatrixDimensions[kCol]);
        }
        else if (iCol >= 0)
        {
            for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
                fResult += m_ppfMatrix[i][iCol] * m_ppfMatrix[i][iCol];
        }
        fResult = sqrt(fResult);
    }
    else 
    {
        if (iRow >= 0)
        {
            for (int j = 0; j < m_aiMatrixDimensions[kRow]; j++)
                fResult += pow(m_ppfMatrix[iRow][j], p);
        }
        else if (iCol >= 0)
        {
            for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
                fResult += pow(m_ppfMatrix[i][iCol], p);
        }
        fResult = pow(fResult,1.F/p);
    }

    return fResult;
}

CMatrix::MatrixError_t CMatrix::setOnes()
{
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::setValue (m_ppfMatrix[i], 1.F, m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

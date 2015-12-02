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

CMatrix::CMatrix( const CMatrix &other ) :
    m_ppfMatrix(0)
{
    assign(other);
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

CMatrix::MatrixError_t CMatrix::transpose()
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
            if (rErr != kMatrixNoError)
                return rErr;
        }
    }

    this->assign(NewMatrix);

    return kMatrixNoError;

}

void CMatrix::assign( const CMatrix &other )
{
    init(other.getNumRows(), other.getNumCols());
    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        other.getRow(i, m_ppfMatrix[i], m_aiMatrixDimensions[kCol]);
}

CMatrix::MatrixError_t CMatrix::mulByElement( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::mulBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

CMatrix::MatrixError_t CMatrix::divByElement( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::divBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

CMatrix::MatrixError_t CMatrix::addByElement( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::addBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

CMatrix::MatrixError_t CMatrix::subByElement( const CMatrix &other )
{
    if (m_aiMatrixDimensions[kRow] != other.getNumRows() || m_aiMatrixDimensions[kCol] != other.getNumCols())
        return kMatrixIllegalFunctionParam;

    for (int i = 0; i < m_aiMatrixDimensions[kRow]; i++)
        CUtil::subBuff(m_ppfMatrix[i], other.getRow(i), m_aiMatrixDimensions[kCol]);

    return kMatrixNoError;
}

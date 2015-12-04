
// standard headers

// project headers
#include "MigEditConfig.h"

#include "ErrorDef.h"

#include "Util.h"
#include "Matrix.h"
#include "NmfIf.h"
#include "Nmf.h"

static const char*  kCMigEditBuildDate             = __DATE__;


//CMigEditIf::CMigEditIf ()
//{
//    // this never hurts
//    this->resetInstance ();
//}
//
//
//CMigEditIf::~CMigEditIf ()
//{
//    this->resetInstance ();
//}

const int  CNmfIf::getVersion (const Version_t eVersionIdx)
{
    int iVersion = 0;

    switch (eVersionIdx)
    {
    case kMajor:
        iVersion    = MigEdit_VERSION_MAJOR; 
        break;
    case kMinor:
        iVersion    = MigEdit_VERSION_MINOR; 
        break;
    case kPatch:
        iVersion    = MigEdit_VERSION_PATCH; 
        break;
    case kNumVersionInts:
        iVersion    = -1;
        break;
    }

    return iVersion;
}
const char*  CNmfIf::getBuildDate ()
{
    return kCMigEditBuildDate;
}

Error_t CNmfIf::create (CNmfIf*& pCMigEdit)
{
    pCMigEdit = new CNmf ();

    if (!pCMigEdit)
        return kUnknownError;


    return kNoError;
}

Error_t CNmfIf::destroy (CNmfIf*& pCMigEdit)
{
    if (!pCMigEdit)
        return kUnknownError;
    
    pCMigEdit->reset ();
    
    delete pCMigEdit;
    pCMigEdit = 0;

    return kNoError;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
CNmfSharedData::CNmfSharedData( int iTemplateLength, int iRankSplit1, int iRankSplit2/*=0*/ ) :
    m_iTemplateLength(iTemplateLength),
    m_fMinError(1e-4F),
    m_iMaxIter(300),
    m_iNumIterations(0),
    m_pfError(0),
    m_bIsInitialized(false)
{
    CUtil::setZero(m_afSparsity, kNumSplits);
    
    for (int i = 0; i < kNumMatrices; i++)
    {
        CUtil::setZero(m_aapCMatrices[i], kNumSplits);
        for (int j = 0; j < kNumSplits; j++)
            m_aabIsUpdated[i][j]    = true;
    }
    m_aiRank[kSplit1]   = iRankSplit1;
    m_aiRank[kSplit2]   = iRankSplit2;

    m_pfError   = new float[m_iMaxIter];
    CUtil::setZero(m_pfError, m_iMaxIter);
}

CNmfSharedData::~CNmfSharedData()
{
    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            delete m_aapCMatrices[i][j];
            m_aapCMatrices[i][j]    = 0;
        }
    }
    if (m_pfError)
        delete [] m_pfError;
    m_pfError   = 0;
}

Error_t CNmfSharedData::setMatrix( Matrices_t eMatrix, MatrixSplit_t eSplit, const CMatrix *pCMatrix )
{
    if (!pCMatrix)
        return kFunctionInvalidArgsError;

    if (eMatrix == kDict)
    {
        if (pCMatrix->getNumCols() != m_aiRank[eSplit])
            return kFunctionInvalidArgsError;
    }
    else if (eMatrix == kAct)
    {
        if (pCMatrix->getNumRows() != m_aiRank[eSplit])
            return kFunctionInvalidArgsError;
    }

    m_aapCMatrices[eMatrix][eSplit] = new CMatrix(*pCMatrix);

    return kNoError;
}

Error_t CNmfSharedData::setIsUpdated( Matrices_t eMatrix, MatrixSplit_t eSplit, bool bIsUpdate /*= true*/ )
{
    m_aabIsUpdated[eMatrix][eSplit] = bIsUpdate;

    return kNoError;
}

Error_t CNmfSharedData::setTerminationCriteria( int iMaxIterations /*= 300*/, float fMinError /*= 1e-4F*/ )
{
    if (iMaxIterations <= 0 || fMinError <= 0)
        return kFunctionInvalidArgsError;

    m_iMaxIter  = iMaxIterations;
    m_fMinError = fMinError;

    return kNoError;
}

Error_t CNmfSharedData::setSparsityLambda( MatrixSplit_t eSplit, float fValue /*= 0*/ )
{
    if (fValue <= 0)
        return kFunctionInvalidArgsError;

    m_afSparsity[eSplit] = fValue;

    return kNoError;
}

Error_t CNmfSharedData::finalizeSettings()
{
    assert(m_aiRank[kSplit1] >= 0 && m_aiRank[kSplit2] >= 0);

    for (int k = 0; k < kNumSplits; k++)
    {
        if (m_aiRank[k] == 0)
        {
            for (int i = 0; i < kNumMatrices; i++)
            {
                if (m_aapCMatrices[i][k])
                    m_aapCMatrices[i][k]->reset();
                m_aabIsUpdated[i][k]    = false;
            }
        }
        else
        {
            if (m_aapCMatrices[kDict][k] != 0)
                if (m_aapCMatrices[kDict][k]->getNumRows() == m_iTemplateLength && m_aapCMatrices[kDict][k]->getNumCols() == m_aiRank[k])
                    continue;
            m_aapCMatrices[kDict][k]->init(m_iTemplateLength, m_aiRank[k]);
            m_aapCMatrices[kDict][k]->setRand();
        }
    }

    if (m_pfError)
        delete [] m_pfError;
    m_pfError   = new float[m_iMaxIter];
    CUtil::setZero(m_pfError, m_iMaxIter);

    m_bIsInitialized    = true;

    return kNoError;
}

Error_t CNmfSharedData::reset()
{
    m_bIsInitialized    = false;
    m_fMinError         = 1e-4F;
    m_iMaxIter          = 300;
    m_iNumIterations    = 0;

    CUtil::setZero(m_afSparsity, kNumSplits);

    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            delete m_aapCMatrices[i][j];
            m_aapCMatrices[i][j]    = 0;
            m_aabIsUpdated[i][j]    = true;
        }
    }

    delete [] m_pfError;
    m_pfError   = 0;

    
    return kNoError;
}

int CNmfSharedData::getMaxIterations() const
{
    return m_iMaxIter;
}

float CNmfSharedData::getMinError() const
{
    return m_fMinError;
}

float CNmfSharedData::getSparsityLambda( MatrixSplit_t eSplit ) const
{
    return m_afSparsity[eSplit];
}

bool CNmfSharedData::getIsUpdated( Matrices_t eMatrix, MatrixSplit_t eSplit ) const
{
    return m_aabIsUpdated[eMatrix][eSplit];
}

CMatrix* CNmfSharedData::getMatrixPtr( Matrices_t eMatrix, MatrixSplit_t eSplit )
{
    return m_aapCMatrices[eMatrix][eSplit];
}

float* CNmfSharedData::getErrorPtr()
{
    return m_pfError;
}

int CNmfSharedData::getResultNumIterations() const
{
    return m_iNumIterations;
}

int CNmfSharedData::getRank( MatrixSplit_t eSplit ) const
{
    return m_aiRank[eSplit];
}

int CNmfSharedData::getTemplateLength() const
{
    return m_iTemplateLength;
}

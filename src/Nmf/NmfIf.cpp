
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
CNmfParametrization::CNmfParametrization(  ) :
    m_iTemplateLength(0),
    m_fMinError(1e-4F),
    m_iMaxIter(300)
{
    reset();
}

CNmfParametrization::~CNmfParametrization()
{
    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            m_aapCMatrices[i][j]    = 0;
        }
    }
}

Error_t CNmfParametrization::setMatrixInit( Matrices_t eMatrix, MatrixSplit_t eSplit, const CMatrix *pCMatrix )
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

Error_t CNmfParametrization::setIsUpdated( Matrices_t eMatrix, MatrixSplit_t eSplit, bool bIsUpdate /*= true*/ )
{
    m_aabIsUpdated[eMatrix][eSplit] = bIsUpdate;

    return kNoError;
}

Error_t CNmfParametrization::setTerminationCriteria( int iMaxIterations /*= 300*/, float fMinError /*= 1e-4F*/ )
{
    if (iMaxIterations <= 0 || fMinError <= 0)
        return kFunctionInvalidArgsError;

    m_iMaxIter  = iMaxIterations;
    m_fMinError = fMinError;

    return kNoError;
}

Error_t CNmfParametrization::setSparsityLambda( MatrixSplit_t eSplit, float fValue /*= 0*/ )
{
    if (fValue <= 0)
        return kFunctionInvalidArgsError;

    m_afSparsity[eSplit] = fValue;

    return kNoError;
}


Error_t CNmfParametrization::reset()
{
    m_fMinError         = 1e-4F;
    m_iMaxIter          = 300;

    CUtil::setZero(m_afSparsity, kNumSplits);
    CUtil::setZero(m_aiRank, kNumSplits);

    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            m_aapCMatrices[i][j]    = 0;
            m_aabIsUpdated[i][j]    = true;
        }
    }
    m_iTemplateLength   = 0;

    return kNoError;
}

int CNmfParametrization::getMaxIterations() const
{
    return m_iMaxIter;
}

float CNmfParametrization::getMinError() const
{
    return m_fMinError;
}

float CNmfParametrization::getSparsityLambda( MatrixSplit_t eSplit ) const
{
    return m_afSparsity[eSplit];
}

bool CNmfParametrization::getIsUpdated( Matrices_t eMatrix, MatrixSplit_t eSplit ) const
{
    return m_aabIsUpdated[eMatrix][eSplit];
}

int CNmfParametrization::getRank( MatrixSplit_t eSplit ) const
{
    return m_aiRank[eSplit];
}

int CNmfParametrization::getTemplateLength() const
{
    return m_iTemplateLength;
}

const CMatrix* CNmfParametrization::getMatrixPtr( Matrices_t eMatrix, MatrixSplit_t eSplit )
{
    return m_aapCMatrices[eMatrix][eSplit];
}

Error_t CNmfParametrization::init( int iTemplateLength, int iRankSplit1, int iRankSplit2 /*= 0*/ )
{
    m_iTemplateLength   = iTemplateLength;
    m_aiRank[kSplit1]   = iRankSplit1;
    m_aiRank[kSplit2]   = iRankSplit2;

    return kNoError;
}

CMatrix* CNmfResult::getMatrixPtr( Matrices_t eMatrix, MatrixSplit_t eSplit )
{
    return m_aapCMatrices[eMatrix][eSplit];
}

CMatrix CNmfResult::getMatrix( Matrices_t eMatrix, MatrixSplit_t eSplit ) const
{
    return *m_aapCMatrices[eMatrix][eSplit];
}

float* CNmfResult::getErrorPtr()
{
    return m_pfError;
}

int CNmfResult::getNumIterations() const
{
    return m_iNumIterations;
}

float CNmfResult::getError() const
{
    if (!m_pfError) 
        return -1.F;

    assert(m_iNumIterations>=0 && m_iNumIterations < m_iMaxIter);

    return m_pfError[m_iNumIterations-1];
}

CNmfResult::CNmfResult() :
    m_iNumIterations(0),
    m_pfError(0),
    m_bIsInitialized(false)
{
    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            m_aapCMatrices[i][j]    = new CMatrix();
        }
    }
    reset ();
}

CNmfResult::~CNmfResult()
{
    reset();
    
    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            delete m_aapCMatrices[i][j];
            m_aapCMatrices[i][j]    = 0;
        }
    }
}

Error_t CNmfResult::init( CNmfParametrization& ParamsAndInit, int iNumObservations )
{
    reset();

    m_iMaxIter  = ParamsAndInit.getMaxIterations();
    m_pfError   = new float[m_iMaxIter];
    CUtil::setZero(m_pfError, m_iMaxIter);
    

    for (int k = 0; k < kNumSplits; k++)
    {
        int iRank = ParamsAndInit.getRank(static_cast<MatrixSplit_t>(k));
        for (int i = 0; i < kNumMatrices; i++)
        {
            if (iRank != 0)
            {
                const CMatrix* pCMat = ParamsAndInit.getMatrixPtr(static_cast<Matrices_t>(i), static_cast<MatrixSplit_t>(k));
                int iNumRows = 0;
                int iNumCols = 0;

                if (i == kDict)
                {
                    iNumRows    = ParamsAndInit.getTemplateLength();
                    iNumCols    = iRank;
                }
                else if (i == kAct)
                {
                    iNumRows    = iRank;
                    iNumCols    = iNumObservations;
                }

                if (!pCMat)
                {
                    m_aapCMatrices[i][k]->reset();
                    m_aapCMatrices[i][k]->init(iNumRows, iNumCols);
                    m_aapCMatrices[i][k]->setRand();
                }
                else
                {
                    if (pCMat->getNumRows() != iNumRows && pCMat->getNumCols() != iNumCols)
                    {
                        m_aapCMatrices[i][k]->reset();
                        m_aapCMatrices[i][k]->init(iNumRows, iNumCols);
                        m_aapCMatrices[i][k]->setRand();
                    }
                    else
                        m_aapCMatrices[i][k] = new CMatrix(*pCMat);
                }                
                if (i == kDict)
                    m_aapCMatrices[i][k]->normalize_I(CMatrix::kPerCol);
            }
        }
    }

    m_bIsInitialized    = true;

    return kNoError;
}

Error_t CNmfResult::reset()
{
    m_bIsInitialized    = false;

    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            m_aapCMatrices[i][j]->reset();
        }
    }
    if (m_pfError)
        delete [] m_pfError;
    m_pfError   = 0;
    m_iMaxIter  = 0;

    return kNoError;
}

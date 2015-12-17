
// standard headers

// project headers
#include "MigEditConfig.h"

#include "ErrorDef.h"

#include "Util.h"
#include "Matrix.h"
#include "NmfIf.h"
#include "Nmf.h"

static const char*  kCMigEditBuildDate             = __DATE__;




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
    m_fNmfMinError(0),
    m_iNmfMaxIter(300),
    m_eMethod(kNoAdaptation),
    m_fAdaptRho(0)
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

Error_t CNmfParametrization::setNmfTerminationCriteria( int iMaxIterations /*= 300*/, float fMinError /*= 1e-3F*/ )
{
    if (iMaxIterations <= 0 || fMinError <= 0)
        return kFunctionInvalidArgsError;

    m_iNmfMaxIter  = iMaxIterations;
    m_fNmfMinError = fMinError;

    return kNoError;
}

Error_t CNmfParametrization::setAdaptTerminationCriteria( int iMaxIterations /*= 300*/, float fMinError /*= 1e-2F*/ )
{
    if (iMaxIterations <= 0 || fMinError <= 0)
        return kFunctionInvalidArgsError;

    m_iAdaptMaxIter  = iMaxIterations;
    m_fAdaptMinError = fMinError;

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
    // default settings
    setNmfTerminationCriteria();
    setAdaptationMethod();
    setAdaptTerminationCriteria();
    setAdaptRho();

    CUtil::setZero(m_afSparsity, kNumSplits);
    CUtil::setZero(m_aiRank, kNumSplits);

    for (int i = 0; i < kNumMatrices; i++)
    {
        for (int j = 0; j < kNumSplits; j++)
        {
            setSparsityLambda((MatrixSplit_t) j);
            m_aapCMatrices[i][j]    = 0;
            m_aabIsUpdated[i][j]    = true;
        }
    }
    m_iTemplateLength   = 0;

    return kNoError;
}

int CNmfParametrization::getNmfMaxIterations() const
{
    return m_iNmfMaxIter;
}

float CNmfParametrization::getNmfMinError() const
{
    return m_fNmfMinError;
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

Error_t CNmfParametrization::setAdaptationMethod( AdaptationMethod_t eMethod  /*= kNoAdaptation*/ )
{
    m_eMethod   = eMethod;

    return kNoError;
}

CNmfParametrization::AdaptationMethod_t CNmfParametrization::getAdaptationMethod() const
{
    return m_eMethod;
}

int CNmfParametrization::getAdaptMaxIterations() const
{
    if (m_eMethod != kNoAdaptation)
        return m_iAdaptMaxIter;
    else
    {
        return 1;
    }
}

float CNmfParametrization::getAdaptMinError() const
{
    return m_fAdaptMinError;
}

float CNmfParametrization::getAdaptRhoThresh() const
{
    return m_fAdaptRho;
}

Error_t CNmfParametrization::setAdaptRho( float fRho /*= .5F*/ )
{
    if (fRho < 0 || fRho > 1.F)
        return kFunctionInvalidArgsError;

    m_fAdaptRho = fRho;

    return kNoError;
}

CNmfResult::CNmfResult() :
    m_iNumNmfIter(0),
    m_iNumAdaptIter(0),
    m_pfNmfError(0),
    m_pfAdaptError(0),
    m_bIsInitialized(false),
    m_iMaxAdaptIter(1),
    m_iMaxNmfIter(0)
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

    m_iMaxNmfIter   = ParamsAndInit.getNmfMaxIterations();
    m_pfNmfError    = new float[m_iMaxNmfIter];
    CUtil::setZero(m_pfNmfError, m_iMaxNmfIter);

    m_iMaxAdaptIter = ParamsAndInit.getAdaptMaxIterations();

    m_pfAdaptError  = new float[m_iMaxAdaptIter];
    CUtil::setZero(m_pfAdaptError, m_iMaxAdaptIter);
 

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
    if (m_pfNmfError)
        delete [] m_pfNmfError;
    m_pfNmfError   = 0;
    m_iMaxNmfIter  = 0;
    if (m_pfAdaptError)
        delete [] m_pfAdaptError;
    m_pfAdaptError   = 0;
    m_iMaxAdaptIter= 0;

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

float* CNmfResult::getNmfErrorPtr()
{
    assert(m_pfNmfError);
    return m_pfNmfError;
}

float* CNmfResult::getAdaptErrorPtr()
{
    assert(m_pfAdaptError);
    return m_pfAdaptError;
}

int CNmfResult::getNumNmfIterations() const
{
    return m_iNumNmfIter;
}

float CNmfResult::getNmfError() const
{
    assert (m_pfNmfError);
    assert(m_iNumNmfIter>=0 && m_iNumNmfIter < m_iMaxNmfIter);

    return m_pfNmfError[m_iNumNmfIter-1];
}

int CNmfResult::getNumAdaptIterations() const
{
    return m_iNumAdaptIter;
}

float CNmfResult::getAdaptError() const
{
    assert (m_pfAdaptError);
    assert(m_iNumAdaptIter>=0 && m_iNumAdaptIter < m_iMaxAdaptIter);

    return m_pfAdaptError[m_iNumNmfIter-1];

}

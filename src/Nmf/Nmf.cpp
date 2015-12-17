
// standard headers
#include <algorithm>

// project headers
#include "ErrorDef.h"

#include "Util.h"
#include "Matrix.h"
#include "Nmf.h"

const float  CNmf::m_kMinOffset = 1e-18F;


CNmf::CNmf () :
    m_phCNmfConfig(0)
{
    // this never hurts
    this->reset ();
}

CNmf::~CNmf ()
{
    this->reset ();
}

Error_t CNmf::init( CNmfParametrization &NmfSharedData )
{
    // set pointers
    m_phCNmfConfig   = &NmfSharedData;

    return kNoError;
}

Error_t CNmf::reset ()
{
    // reset buffers and variables to default values
    m_phCNmfConfig   = 0;

    return kNoError;
}

Error_t CNmf::process( const CMatrix *pCInput, CNmfResult& NmfResult )
{
    if (!pCInput)
        return kFunctionInvalidArgsError;

    // init and configuration stuff
    int aiRank[kNumSplits]          = {0,0};
    int iMaxIter                    = m_phCNmfConfig->getAdaptMaxIterations();
    CMatrix *aaCMatrix[kNumMatrices][kNumSplits];
    
    // result stuff
    float *phfErr                   = 0;

    NmfResult.init(*m_phCNmfConfig, pCInput->getNumCols());
    phfErr  = NmfResult.getAdaptErrorPtr();

    for (int j = 0; j < kNumSplits; j++)
        aiRank[j]   = m_phCNmfConfig->getRank((MatrixSplit_t)j);

    for (int i = 0; i < kNumMatrices; i++)
        for (int j = 0; j < kNumSplits; j++)
        {
            aaCMatrix[i][j]  = NmfResult.getMatrixPtr((Matrices_t)i, (MatrixSplit_t)j);
            if (!aaCMatrix[i][j])
                return kUnknownError;
        }

    if (m_phCNmfConfig->getAdaptationMethod() == CNmfParametrization::kNoAdaptation)
        iMaxIter    = m_phCNmfConfig->getAdaptMaxIterations();

    // verify that not both ranks are 0
    if (aiRank[kSplit1] <= 0 && aiRank[kSplit2] <= 0)
        return kFunctionInvalidArgsError;
    
    // at least one pair of template and activation must exist
    if ((!aaCMatrix[kDict][kSplit2] && !aaCMatrix[kDict][kSplit1])  || (!aaCMatrix[kAct][kSplit2] && !aaCMatrix[kAct][kSplit1]))
        return kFunctionInvalidArgsError;

    // verify that matrix dimensions make sense
    if ((aaCMatrix[kDict][kSplit1]->getNumRows()+aaCMatrix[kDict][kSplit2]->getNumRows() != pCInput->getNumRows())  || 
        (std::max(aaCMatrix[kAct][kSplit1]->getNumCols(),aaCMatrix[kAct][kSplit2]->getNumCols()) != pCInput->getNumCols()))
        return kFunctionInvalidArgsError;

    for (int k = 0; k < iMaxIter; k++)
    {
        if (m_phCNmfConfig->getAdaptationMethod() != CNmfParametrization::kNoAdaptation)
        {
            //assume drums are in split1
            m_phCNmfConfig->setIsUpdated(kDict, kSplit1, false);

            m_phCNmfConfig->setIsUpdated(kDict, kSplit2, true);
            NmfResult.getMatrixPtr(kDict, kSplit2)->setRand();
            m_phCNmfConfig->setIsUpdated(kAct, kSplit1, true);
            NmfResult.getMatrixPtr(kAct, kSplit1)->setRand();
            m_phCNmfConfig->setIsUpdated(kAct, kSplit2, true);
            NmfResult.getMatrixPtr(kAct, kSplit2)->setRand();
        }
        runNmf(pCInput, NmfResult);

        phfErr[k]   = NmfResult.getNmfError();

        if (isRelativeErrorBelowThresh(k, phfErr))
        {
            NmfResult.m_iNumAdaptIter   = k+1;
            break;
        }

        if (m_phCNmfConfig->getAdaptationMethod() == CNmfParametrization::kAdaptation01)
        {
            float fAdaptCoeff = 1.F/(1<<k);
            adaptTemplateXCorr(aaCMatrix, aiRank, fAdaptCoeff);

            m_phCNmfConfig->setIsUpdated(kDict, kSplit1, false);
            m_phCNmfConfig->setIsUpdated(kDict, kSplit2, false);

            m_phCNmfConfig->setIsUpdated(kAct, kSplit1, true);
            NmfResult.getMatrixPtr(kAct, kSplit1)->setRand();
            m_phCNmfConfig->setIsUpdated(kAct, kSplit2, true);
            NmfResult.getMatrixPtr(kAct, kSplit2)->setRand();

            runNmf(pCInput, NmfResult);

            phfErr[k]   = NmfResult.getNmfError();

            if (isRelativeErrorBelowThresh(k, phfErr))
            {
                NmfResult.m_iNumAdaptIter   = k+1;
                break;
            }
        }
        else if (m_phCNmfConfig->getAdaptationMethod() == CNmfParametrization::kAdaptation02)
        {
            //assume drums are in split1
            m_phCNmfConfig->setIsUpdated(kDict, kSplit1, false);
            m_phCNmfConfig->setIsUpdated(kAct, kSplit1, false);

            m_phCNmfConfig->setIsUpdated(kDict, kSplit2, true);
            NmfResult.getMatrixPtr(kDict, kSplit2)->setRand();
            m_phCNmfConfig->setIsUpdated(kAct, kSplit2, true);
            NmfResult.getMatrixPtr(kAct, kSplit2)->setRand();
        }
    }

    return kNoError;
}

bool CNmf::isRelativeErrorBelowThresh( int iCurrIteration, const float *pfErr ) const
{
    assert (pfErr);

    if (iCurrIteration < 1)
        return false;

    return (abs(pfErr[iCurrIteration] - pfErr[iCurrIteration-1])/(pfErr[0] - pfErr[iCurrIteration] + m_kMinOffset)) < m_phCNmfConfig->getNmfMinError();
}

void CNmf::runNmf( const CMatrix * pCInput, CNmfResult &NmfResult )
{
    CMatrix CXHat;
    CMatrix *aaCMatrix[kNumMatrices][kNumSplits];
    int aiRank[kNumSplits]          = {m_phCNmfConfig->getRank(kSplit1), m_phCNmfConfig->getRank(kSplit2)};
    float afWeight[kNumSplits]      = {1.F,1.F};
    float afSparsity[kNumSplits]    = {m_phCNmfConfig->getSparsityLambda(kSplit1), m_phCNmfConfig->getSparsityLambda(kSplit2)};
    float *phfErr                   = NmfResult.getNmfErrorPtr();
    int iMaxIter                    = m_phCNmfConfig->getNmfMaxIterations();

    for (int i = 0; i < kNumMatrices; i++)
        for (int j = 0; j < kNumSplits; j++)
            aaCMatrix[i][j]  = NmfResult.getMatrixPtr((Matrices_t)i, (MatrixSplit_t)j);

    // compute weighting to weight the small split higher
    if (aiRank[kSplit1]*aiRank[kSplit2] != 0)
    {
        if (aiRank[kSplit1] > aiRank[kSplit2])
        {
            afWeight[kSplit2] = (aiRank[kSplit2]  + aiRank[kSplit1]) * 1.F / aiRank[kSplit2]; //alpha
            afWeight[kSplit1] = aiRank[kSplit1] * 1.F/ (aiRank[kSplit2]  + aiRank[kSplit1]); //beta
        }
        else
        {
            afWeight[kSplit1] = (aiRank[kSplit1]  + aiRank[kSplit2]) * 1.F / aiRank[kSplit1]; //alpha
            afWeight[kSplit2] = aiRank[kSplit2] * 1.F/ (aiRank[kSplit1]  + aiRank[kSplit2]); //beta
        }
        CXHat   = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1] * afWeight[kSplit1] + *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2] * afWeight[kSplit2];
    }
    else if (aiRank[kSplit2] > 0)
        CXHat   = *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2];
    else if (aiRank[kSplit1] > 0)
        CXHat   = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1];

    CXHat.addC_I(m_kMinOffset);

    for (int k = 0; k < iMaxIter; k++)
    {
        float fSparsityCost = 0;

        // (X./approx)
        CXHat = pCInput->divByElement(CXHat);

        for (int j = 0; j < kNumSplits; j++)
        {
            if (aiRank[j] <= 0)
                continue;
            // update rules
            // HH = HH .* ((beta * WH)'* (X./approx))./((beta * WH)'*rep);
            if (m_phCNmfConfig->getIsUpdated(kAct, static_cast<MatrixSplit_t>(j)))            
            {
                aaCMatrix[kAct][j]->mulByElement_I(
                    (aaCMatrix[kDict][j]->transpose().mulC_I(afWeight[j]) * CXHat ).divByElement_I(
                    aaCMatrix[kDict][j]->transpose().mulC_I(afWeight[j]).mulByOnes(pCInput->getNumRows(),pCInput->getNumCols()).addC_I(afSparsity[j])));
                aaCMatrix[kAct][j]->setZeroBelowThresh(m_kMinOffset*aaCMatrix[kAct][kSplit1]->getMax());
            }

            if (m_phCNmfConfig->getIsUpdated(kDict, static_cast<MatrixSplit_t>(j)))            
            {
                // WH = WH .* ((X./approx)*(beta * HH)')./(rep*(beta * HH)');
                aaCMatrix[kDict][j]->mulByElement_I(
                    (CXHat * aaCMatrix[kAct][j]->transpose().mulC_I(afWeight[j])).divByElement_I(
                    ((aaCMatrix[kAct][j]->mulC_I(afWeight[j])).mulByOnes(pCInput->getNumCols(),pCInput->getNumRows())).transpose()));

                // normalization
                //aaCMatrix[kDict][kSplit1]->setZeroBelowThresh(m_kMinOffset);
                aaCMatrix[kDict][j]->normalize_I(CMatrix::kPerCol);
                aaCMatrix[kDict][j]->setZeroBelowThresh(m_kMinOffset*aaCMatrix[kAct][kSplit1]->getMax());
            }
        }
 
        if (aiRank[kSplit1] > 0 && aiRank[kSplit2] > 0)
            CXHat           = (*aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1]).mulC_I(afWeight[kSplit1]).addByElement_I(
            (*aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2]).mulC_I(afWeight[kSplit2]));
        else if (aiRank[kSplit2] > 0)
            CXHat           = *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2];
        else if (aiRank[kSplit1] > 0)
            CXHat           = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1];
        CXHat.addC_I(m_kMinOffset);

        for (int j = 0; j < kNumSplits; j++)
            fSparsityCost  += afSparsity[j] * aaCMatrix[kAct][j]->getNorm();

        if (fSparsityCost > 0)
            phfErr[k]  = CMatrix::getKlDivergence(*pCInput, CXHat + fSparsityCost);
        else 
            phfErr[k]  = CMatrix::getKlDivergence(*pCInput, CXHat);

        // termination criteria
        if (isRelativeErrorBelowThresh(k, phfErr))
        {
            for (int i = 0; i < kNumSplits; i++)
            {
                if (aiRank[i] > 0)
                    aaCMatrix[kAct][kSplit1]->setZeroBelowThresh(m_kMinOffset*aaCMatrix[kAct][kSplit1]->getMax());
            }

            NmfResult.m_iNumNmfIter    = k+1;
            break;
        }
    }
}

void CNmf::adaptTemplateXCorr( CMatrix * paaCMatrix[kNumMatrices][kNumSplits], const int aiRank[], float fAdaptCoeff )
{
    //normalized cross correlation
    CMatrix Rho = calcCrossCorr(paaCMatrix[kAct]);
        
    for (int i = 0; i < aiRank[kSplit1]; i++)
    {
        int iNumRows            = paaCMatrix[kDict][kSplit1]->getNumRows();
        int iNumCorrelatedElems = 0;
        float *pfAdapt          = new float [iNumRows];
        CUtil::setZero(pfAdapt, iNumRows);
        for (int j = 0; j < aiRank[kSplit2]; j++)
        {
            if (Rho.getElement(i,j) >= m_phCNmfConfig->getAdaptRhoThresh())
            {
                for (int k = 0; k < iNumRows; k++)
                {
                    pfAdapt[k] += fAdaptCoeff * paaCMatrix[kDict][kSplit2]->getElement(k,j) * Rho.getElement(i,j);
                    paaCMatrix[kDict][kSplit2]->setElement(k,j,std::rand()*1.F/RAND_MAX);
                }
                iNumCorrelatedElems++;
            }
        }
        CUtil::mulBuffC(pfAdapt, 1.F/iNumCorrelatedElems, iNumRows);

        for (int k = 0; k < iNumRows; k++)
            paaCMatrix[kDict][kSplit1]->setElement(k,i, (1-fAdaptCoeff) * paaCMatrix[kDict][kSplit1]->getElement(k,i) + pfAdapt[k]);
    }
}

CMatrix CNmf::calcCrossCorr( CMatrix * paCMatrix[kNumSplits] )
{
    int iNumRows = paCMatrix[kSplit1]->getNumRows();
    int iNumCols = paCMatrix[kSplit2]->getNumRows();

    CMatrix Rho(iNumRows, iNumCols);
  
    for (int i = 0; i < iNumRows; i++)
    {
        for (int j = 0; j < iNumCols; j++)
        {
            float fRho = CUtil::mulBuffScalar(paCMatrix[kSplit1]->getRowPtr(i), paCMatrix[kSplit2]->getRowPtr(j), paCMatrix[kSplit1]->getNumCols());
            fRho   /= (paCMatrix[kSplit1]->getRowNorm(i,2) * paCMatrix[kSplit2]->getRowNorm(j,2));
            Rho.setElement(i,j, fRho);
        }
    }
    
    return Rho;
}


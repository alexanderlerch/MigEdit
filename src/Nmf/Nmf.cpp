
// standard headers

// project headers
#include "ErrorDef.h"

#include "Util.h"
#include "Matrix.h"
#include "Nmf.h"

const float  CNmf::m_kMinOffset = 1e-18F;


CNmf::CNmf () :
    m_phfErr(0),
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
    m_phfErr                = 0;
    m_phCNmfConfig   = 0;

    return kNoError;
}

Error_t CNmf::process( const CMatrix *pCInput, CNmfResult& NmfResult )
{
    if (!pCInput)
        return kFunctionInvalidArgsError;

    CMatrix CXHat;
    CMatrix *aaCMatrix[kNumMatrices][kNumSplits];
    int aiRank[kNumSplits]          = {0,0};
    float afWeight[kNumSplits]      = {1.F,1.F};
    float afSparsity[kNumSplits]    = {0,0};
    int iMaxIter                    = m_phCNmfConfig->getMaxIterations();


    NmfResult.init(*m_phCNmfConfig, pCInput->getNumCols());
    m_phfErr                        = NmfResult.getErrorPtr();

    for (int i = 0; i < kNumMatrices; i++)
        for (int j = 0; j < kNumSplits; j++)
            aaCMatrix[i][j]  = NmfResult.getMatrixPtr((Matrices_t)i, (MatrixSplit_t)j);

    for (int j = 0; j < kNumSplits; j++)
    {
        afSparsity[j]   = m_phCNmfConfig->getSparsityLambda((MatrixSplit_t)j);
        
        aiRank[j]       = m_phCNmfConfig->getRank((MatrixSplit_t)j);
    }
    
    if ((!aaCMatrix[kDict][kSplit2] && !aaCMatrix[kDict][kSplit1])  || (!aaCMatrix[kAct][kSplit2] && !aaCMatrix[kAct][kSplit1]))
        return kFunctionInvalidArgsError;

    if (aiRank[kSplit1] <= 0 && aiRank[kSplit2] <= 0)
        return kFunctionInvalidArgsError;
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
    }

    if (aiRank[kSplit1] > 0 && aiRank[kSplit2] > 0)
        CXHat           = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1] * afWeight[kSplit1] + *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2] * afWeight[kSplit2] + m_kMinOffset;
    else if (aiRank[kSplit2] > 0)
        CXHat           = *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2] * afWeight[kSplit2] + m_kMinOffset;
    else if (aiRank[kSplit1] > 0)
        CXHat           = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1] * afWeight[kSplit1] + m_kMinOffset;

    for (int k = 0; k < iMaxIter; k++)
    {
        float fSparsityCost = 0;

        CXHat = pCInput->divByElement(CXHat);

        if (aiRank[kSplit2] > 0)
        {
            // update rules
            // HD = HD .* ((alpha * WD)'* (X./approx))./((alpha * WD)'*rep + sparsity);
            if (m_phCNmfConfig->getIsUpdated(kAct, kSplit2))            
                aaCMatrix[kAct][kSplit2]->mulByElement_I((aaCMatrix[kDict][kSplit2]->transpose() * CXHat * afWeight[kSplit2]).divByElement(aaCMatrix[kDict][kSplit2]->transpose().mulByOnes(pCInput->getNumRows(),pCInput->getNumCols())*afWeight[kSplit2] + afSparsity[kSplit2]));

            // WD = WD .* ((X./approx)*(alpha * HD)')./(rep*(alpha * HD)');
            if (m_phCNmfConfig->getIsUpdated(kDict, kSplit2))            
            {
                aaCMatrix[kDict][kSplit2]->mulByElement_I((CXHat * aaCMatrix[kAct][kSplit2]->transpose() * afWeight[kSplit2]).divByElement((aaCMatrix[kAct][kSplit2]->mulByOnes(pCInput->getNumCols(),pCInput->getNumRows()) * afWeight[kSplit2]).transpose()));
                //aaCMatrix[kDict][kSplit2]->setZeroBelowThresh(m_kMinOffset);
                // normalization
                aaCMatrix[kDict][kSplit2]->normalize_I(CMatrix::kPerCol);
            }
        }
        if (aiRank[kSplit1] > 0)
        {
            // update rules
            // HH = HH .* ((beta * WH)'* (X./approx))./((beta * WH)'*rep);
            if (m_phCNmfConfig->getIsUpdated(kAct, kSplit2))            
                aaCMatrix[kAct][kSplit1]->mulByElement_I((aaCMatrix[kDict][kSplit1]->transpose() * CXHat * afWeight[kSplit1]).divByElement(aaCMatrix[kDict][kSplit1]->transpose().mulByOnes(pCInput->getNumRows(),pCInput->getNumCols())*afWeight[kSplit1] + afSparsity[kSplit1]));

            if (m_phCNmfConfig->getIsUpdated(kDict, kSplit1))            
            {
                // WH = WH .* ((X./approx)*(beta * HH)')./(rep*(beta * HH)');
                aaCMatrix[kDict][kSplit1]->mulByElement_I((CXHat * aaCMatrix[kAct][kSplit1]->transpose() * afWeight[kSplit1]).divByElement((aaCMatrix[kAct][kSplit1]->mulByOnes(pCInput->getNumCols(),pCInput->getNumRows()) * afWeight[kSplit1]).transpose()));
                //aaCMatrix[kDict][kSplit1]->setZeroBelowThresh(m_kMinOffset);
                // normalization
                aaCMatrix[kDict][kSplit1]->normalize_I(CMatrix::kPerCol);
            }
        }
        if (aiRank[kSplit1] > 0 && aiRank[kSplit2] > 0)
        {
            CXHat           = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1] * afWeight[kSplit1] + *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2] * afWeight[kSplit2] + m_kMinOffset;
            fSparsityCost   = afSparsity[kSplit2] * aaCMatrix[kAct][kSplit2]->getSum() + afSparsity[kSplit1] * aaCMatrix[kAct][kSplit1]->getSum();
        }
        else if (aiRank[kSplit2] > 0)
        {
            CXHat           = *aaCMatrix[kDict][kSplit2] * *aaCMatrix[kAct][kSplit2] * afWeight[kSplit2] + m_kMinOffset;
            fSparsityCost   = afSparsity[kSplit2] * aaCMatrix[kAct][kSplit2]->getSum();
        }
        else if (aiRank[kSplit1] > 0)
        {
            CXHat           = *aaCMatrix[kDict][kSplit1] * *aaCMatrix[kAct][kSplit1] * afWeight[kSplit1] + m_kMinOffset;
            fSparsityCost   = afSparsity[kSplit1] * aaCMatrix[kAct][kSplit1]->getSum();
        }


        if (fSparsityCost > 0)
            m_phfErr[k]  = calcKlDivergence(*pCInput, CXHat + fSparsityCost);
        else 
            m_phfErr[k]  = calcKlDivergence(*pCInput, CXHat);

        // termination criteria
        if (isRelativeErrorBelowThresh(k))
        {
            for (int i = 0; i < kNumSplits; i++)
            {
                if (aiRank[i] > 0)
                    aaCMatrix[kAct][kSplit1]->setZeroBelowThresh(m_kMinOffset*aaCMatrix[kAct][kSplit1]->getMax());
            }
            
            NmfResult.m_iNumIterations    = k;
            break;
        }
    }

    return kNoError;
}

bool CNmf::isRelativeErrorBelowThresh( int iCurrIteration ) const
{
    if (iCurrIteration < 1)
        return false;
    return (abs(m_phfErr[iCurrIteration] - m_phfErr[iCurrIteration-1])/(m_phfErr[0] - m_phfErr[iCurrIteration] + m_kMinOffset) < m_phCNmfConfig->getMinError());
}

float CNmf::calcKlDivergence( const CMatrix &Mat1, const CMatrix &Mat2 ) const
{
    int iNumRows    = Mat1.getNumRows();
    int iNumCols    = Mat1.getNumCols();
    float fResult   = 0;
    assert (iNumRows == Mat2.getNumRows() || iNumCols == Mat2.getNumCols());
    
    for (int i = 0; i < iNumRows; i++)
    {
        for (int j = 0; j < iNumCols; j++)
        {
            float fElem1 = Mat1.getElement(i,j);
            float fElem2 = Mat2.getElement(i,j);
            fResult     += fElem1 * (std::log(fElem1 + m_kMinOffset) - std::log(fElem2 + m_kMinOffset)) - fElem1 + fElem2;
        }
    }

    return fResult;
}


// standard headers

// project headers
#include "ErrorDef.h"

#include "Util.h"
#include "Matrix.h"
#include "Nmf.h"

const float  CNmf::m_kMinOffset = 1e-18F;


CNmf::CNmf () :
    m_pfErr(0),
    m_pCConfigAndResults(0)
{
    // this never hurts
    this->reset ();
}

CNmf::~CNmf ()
{
    this->reset ();
}

Error_t CNmf::init( CNmfSharedData &NmfSharedData )
{
    // set pointers
    m_pfErr = NmfSharedData.getErrorPtr();
    m_pCConfigAndResults  = &NmfSharedData;

    return kNoError;
}

Error_t CNmf::reset ()
{
    // reset buffers and variables to default values
    m_pfErr         = 0;
    m_pCConfigAndResults  = 0;

    return kNoError;
}

Error_t CNmf::process( const CMatrix *pCInput )
{
    if (!pCInput)
        return kFunctionInvalidArgsError;

    CMatrix CXHat;
    CMatrix *aaMatrix[kNumMatrices][kNumSplits];
    int aiRank[kNumSplits]          = {0,0};
    float m_afWeight[kNumSplits]    = {1.F,1.F};
    float afSparsity[kNumSplits]    = {0,0};
    int iMaxIter                    = m_pCConfigAndResults->getMaxIterations();

    for (int i = 0; i < kNumMatrices; i++)
        for (int j = 0; j < kNumSplits; j++)
            aaMatrix[i][j]  = m_pCConfigAndResults->getMatrixPtr((Matrices_t)i, (MatrixSplit_t)j);

    for (int j = 0; j < kNumSplits; j++)
    {
        afSparsity[j]   = m_pCConfigAndResults->getSparsityLambda((MatrixSplit_t)j);
        
        aiRank[j]       = m_pCConfigAndResults->getRank((MatrixSplit_t)j);
        if (aiRank[j])
        {
            if (aaMatrix[kDict][j])
            {
                if (pCInput->getNumRows() != aaMatrix[kDict][j]->getNumRows())
                {
                    aaMatrix[kDict][j]->init(m_pCConfigAndResults->getTemplateLength(), aiRank[j]);
                    aaMatrix[kDict][j]->setRand();
                }
            }
            else
            {
                aaMatrix[kDict][j]->init(m_pCConfigAndResults->getTemplateLength(), aiRank[j]);
                aaMatrix[kDict][j]->setRand();
            }

            // normalization of template matrices
            aaMatrix[kDict][j]->normalize_I(CMatrix::kPerCol);

            if (aaMatrix[kAct][j])
            {
                if (pCInput->getNumCols() != aaMatrix[kAct][j]->getNumCols())
                {
                    aaMatrix[kAct][j]->init(aaMatrix[kDict][j]->getNumCols(), pCInput->getNumCols());
                    aaMatrix[kAct][j]->setRand();
                }
            }
            else
            {
                aaMatrix[kAct][j]->init(aaMatrix[kDict][j]->getNumCols(), pCInput->getNumCols());
                aaMatrix[kAct][j]->setRand();
            }
        }
    }
    if (aiRank[kSplit1] <= 0 && aiRank[kSplit2] <= 0)
        return kFunctionInvalidArgsError;
    if ((!aaMatrix[kDict][kSplit2] && !aaMatrix[kDict][kSplit1])  || (!aaMatrix[kAct][kSplit2] && !aaMatrix[kAct][kSplit1]))
        return kFunctionInvalidArgsError;

    if (aiRank[kSplit1]*aiRank[kSplit2] != 0)
    {
        if (aiRank[kSplit1] > aiRank[kSplit2])
        {
            m_afWeight[kSplit2] = (aiRank[kSplit2]  + aiRank[kSplit1]) * 1.F / aiRank[kSplit2]; //alpha
            m_afWeight[kSplit1] = aiRank[kSplit1] * 1.F/ (aiRank[kSplit2]  + aiRank[kSplit1]); //beta
        }
        else
        {
            m_afWeight[kSplit1] = (aiRank[kSplit1]  + aiRank[kSplit2]) * 1.F / aiRank[kSplit1]; //alpha
            m_afWeight[kSplit2] = aiRank[kSplit2] * 1.F/ (aiRank[kSplit1]  + aiRank[kSplit2]); //beta
        }
    }

    if (aiRank[kSplit1] > 0 && aiRank[kSplit2] > 0)
        CXHat           = *aaMatrix[kDict][kSplit1] * *aaMatrix[kAct][kSplit1] * m_afWeight[kSplit1] + *aaMatrix[kDict][kSplit2] * *aaMatrix[kAct][kSplit2] * m_afWeight[kSplit2] + m_kMinOffset;
    else if (aiRank[kSplit2] > 0)
        CXHat           = *aaMatrix[kDict][kSplit2] * *aaMatrix[kAct][kSplit2] * m_afWeight[kSplit2] + m_kMinOffset;
    else if (aiRank[kSplit1] > 0)
        CXHat           = *aaMatrix[kDict][kSplit1] * *aaMatrix[kAct][kSplit1] * m_afWeight[kSplit1] + m_kMinOffset;

    for (int k = 0; k < iMaxIter; k++)
    {
        float fSparsityCost = 0;

        CXHat = pCInput->divByElement(CXHat);

        if (aiRank[kSplit2] > 0)
        {
            // update rules
            // WD = WD .* ((X./approx)*(alpha * HD)')./(rep*(alpha * HD)');
            if (m_pCConfigAndResults->getIsUpdated(kDict, kSplit2))            
                aaMatrix[kDict][kSplit2]->mulByElement_I((CXHat * aaMatrix[kAct][kSplit2]->transpose() * m_afWeight[kSplit2]).divByElement((aaMatrix[kAct][kSplit2]->mulByOnes(pCInput->getNumCols(),pCInput->getNumRows()) * m_afWeight[kSplit2]).transpose()));
            // HD = HD .* ((alpha * WD)'* (X./approx))./((alpha * WD)'*rep + sparsity);
            if (m_pCConfigAndResults->getIsUpdated(kAct, kSplit2))            
                aaMatrix[kAct][kSplit2]->mulByElement_I((aaMatrix[kDict][kSplit2]->transpose() * CXHat * m_afWeight[kSplit2]).divByElement(aaMatrix[kDict][kSplit2]->transpose().mulByOnes(pCInput->getNumRows(),pCInput->getNumCols())*m_afWeight[kSplit2] + afSparsity[kSplit2]));

            aaMatrix[kDict][kSplit2]->normalize_I(CMatrix::kPerCol);
        }
        if (aiRank[kSplit1] > 0)
        {
            // update rules
            // WH = WH .* ((X./approx)*(beta * HH)')./(rep*(beta * HH)');
            if (m_pCConfigAndResults->getIsUpdated(kDict, kSplit1))            
                aaMatrix[kDict][kSplit1]->mulByElement_I((CXHat * aaMatrix[kAct][kSplit1]->transpose() * m_afWeight[kSplit1]).divByElement((aaMatrix[kAct][kSplit1]->mulByOnes(pCInput->getNumCols(),pCInput->getNumRows()) * m_afWeight[kSplit1]).transpose()));
            // HH = HH .* ((beta * WH)'* (X./approx))./((beta * WH)'*rep);
            if (m_pCConfigAndResults->getIsUpdated(kAct, kSplit2))            
                aaMatrix[kAct][kSplit1]->mulByElement_I((aaMatrix[kDict][kSplit1]->transpose() * CXHat * m_afWeight[kSplit1]).divByElement(aaMatrix[kDict][kSplit1]->transpose().mulByOnes(pCInput->getNumRows(),pCInput->getNumCols())*m_afWeight[kSplit1] + afSparsity[kSplit1]));

            // normalization
            aaMatrix[kDict][kSplit1]->normalize_I(CMatrix::kPerCol);
        }
        if (aiRank[kSplit1] > 0 && aiRank[kSplit2] > 0)
        {
            CXHat           = *aaMatrix[kDict][kSplit1] * *aaMatrix[kAct][kSplit1] * m_afWeight[kSplit1] + *aaMatrix[kDict][kSplit2] * *aaMatrix[kAct][kSplit2] * m_afWeight[kSplit2] + m_kMinOffset;
            fSparsityCost   = afSparsity[kSplit2] * aaMatrix[kAct][kSplit2]->getSum() + afSparsity[kSplit1] * aaMatrix[kAct][kSplit1]->getSum();
        }
        else if (aiRank[kSplit2] > 0)
        {
            CXHat           = *aaMatrix[kDict][kSplit2] * *aaMatrix[kAct][kSplit2] * m_afWeight[kSplit2] + m_kMinOffset;
            fSparsityCost   = afSparsity[kSplit2] * aaMatrix[kAct][kSplit2]->getSum();
        }
        else if (aiRank[kSplit1] > 0)
        {
            CXHat           = *aaMatrix[kDict][kSplit1] * *aaMatrix[kAct][kSplit1] * m_afWeight[kSplit1] + m_kMinOffset;
            fSparsityCost   = afSparsity[kSplit1] * aaMatrix[kAct][kSplit1]->getSum();
        }


        if (fSparsityCost > 0)
            m_pfErr[k]  = calcKlDivergence(*pCInput, CXHat + fSparsityCost);
        else 
            m_pfErr[k]  = calcKlDivergence(*pCInput, CXHat);

        // termination criteria
        if (isRelativeErrorBelowThresh(k))
        {
            m_pCConfigAndResults->m_iNumIterations    = k;
            break;
        }
    }

    return kNoError;
}

bool CNmf::isRelativeErrorBelowThresh( int iCurrIteration ) const
{
    if (iCurrIteration < 1)
        return false;
    return (abs(m_pfErr[iCurrIteration] - m_pfErr[iCurrIteration-1])/(m_pfErr[0] - m_pfErr[iCurrIteration] + m_kMinOffset) < m_pCConfigAndResults->getMinError());
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

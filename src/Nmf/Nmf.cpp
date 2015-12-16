
// standard headers
#include <algorithm>

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
    
    // at least one pair of template and activation must exist
    if ((!aaCMatrix[kDict][kSplit2] && !aaCMatrix[kDict][kSplit1])  || (!aaCMatrix[kAct][kSplit2] && !aaCMatrix[kAct][kSplit1]))
        return kFunctionInvalidArgsError;

    // verify that matrix dimensions make sense
    if ((aaCMatrix[kDict][kSplit1]->getNumRows()+aaCMatrix[kDict][kSplit2]->getNumRows() != pCInput->getNumRows())  || 
        (std::max(aaCMatrix[kAct][kSplit1]->getNumCols(),aaCMatrix[kAct][kSplit2]->getNumCols()) != pCInput->getNumCols()))
        return kFunctionInvalidArgsError;

    // verify that not both ranks are 0
    if (aiRank[kSplit1] <= 0 && aiRank[kSplit2] <= 0)
        return kFunctionInvalidArgsError;

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
                    aaCMatrix[kAct][j]->mulByElement_I(
                    (aaCMatrix[kDict][j]->transpose().mulC_I(afWeight[j]) * CXHat ).divByElement_I(
                    aaCMatrix[kDict][j]->transpose().mulC_I(afWeight[j]).mulByOnes(pCInput->getNumRows(),pCInput->getNumCols()).addC_I(afSparsity[j])));

            if (m_phCNmfConfig->getIsUpdated(kDict, static_cast<MatrixSplit_t>(j)))            
            {
                // WH = WH .* ((X./approx)*(beta * HH)')./(rep*(beta * HH)');
                aaCMatrix[kDict][j]->mulByElement_I(
                    (CXHat * aaCMatrix[kAct][j]->transpose().mulC_I(afWeight[j])).divByElement_I(
                    ((aaCMatrix[kAct][j]->mulC_I(afWeight[j])).mulByOnes(pCInput->getNumCols(),pCInput->getNumRows())).transpose()));

                // normalization
                //aaCMatrix[kDict][kSplit1]->setZeroBelowThresh(m_kMinOffset);
                aaCMatrix[kDict][j]->normalize_I(CMatrix::kPerCol);
            }
        }
        //if (aiRank[kSplit2] > 0)
        //{
        //    // update rules
        //    // HD = HD .* ((alpha * WD)'* (X./approx))./((alpha * WD)'*rep + sparsity);
        //    if (m_phCNmfConfig->getIsUpdated(kAct, kSplit2))            
        //        aaCMatrix[kAct][kSplit2]->mulByElement_I((aaCMatrix[kDict][kSplit2]->transpose() * CXHat * afWeight[kSplit2]).divByElement(aaCMatrix[kDict][kSplit2]->transpose().mulByOnes(pCInput->getNumRows(),pCInput->getNumCols())*afWeight[kSplit2] + afSparsity[kSplit2]));

        //    // WD = WD .* ((X./approx)*(alpha * HD)')./(rep*(alpha * HD)');
        //    if (m_phCNmfConfig->getIsUpdated(kDict, kSplit2))            
        //    {
        //        aaCMatrix[kDict][kSplit2]->mulByElement_I((CXHat * aaCMatrix[kAct][kSplit2]->transpose() * afWeight[kSplit2]).divByElement((aaCMatrix[kAct][kSplit2]->mulByOnes(pCInput->getNumCols(),pCInput->getNumRows()) * afWeight[kSplit2]).transpose()));
        //        //aaCMatrix[kDict][kSplit2]->setZeroBelowThresh(m_kMinOffset);
        //        // normalization
        //        aaCMatrix[kDict][kSplit2]->normalize_I(CMatrix::kPerCol);
        //    }
        //}
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
            m_phfErr[k]  = CMatrix::getKlDivergence(*pCInput, CXHat + fSparsityCost);
        else 
            m_phfErr[k]  = CMatrix::getKlDivergence(*pCInput, CXHat);

        // termination criteria
        if (isRelativeErrorBelowThresh(k))
        {
            for (int i = 0; i < kNumSplits; i++)
            {
                if (aiRank[i] > 0)
                    aaCMatrix[kAct][kSplit1]->setZeroBelowThresh(m_kMinOffset*aaCMatrix[kAct][kSplit1]->getMax());
            }
            
            NmfResult.m_iNumIterations    = k+1;
            break;
        }
    }

    return kNoError;
}

bool CNmf::isRelativeErrorBelowThresh( int iCurrIteration ) const
{
    if (iCurrIteration < 1)
        return false;
    return (abs(m_phfErr[iCurrIteration] - m_phfErr[iCurrIteration-1])/(m_phfErr[0] - m_phfErr[iCurrIteration] + m_kMinOffset)) < m_phCNmfConfig->getMinError();
}


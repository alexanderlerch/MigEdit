#include <algorithm>
#include <cstdlib>

#include "Util.h"
#include "AudioInfo.h"
#include "Preproc.h"

CPreproc::CPreproc() :
    m_phAudioInfo(0),
    m_bIsInitialized(false)
{
    for (int i = 0; i < kNumPreprocSteps; i++)
        m_abIsStepActive[i] = false;

    reset ();
}

Error_t CPreproc::create( CPreproc*& pCPreproc )
{
    pCPreproc = new CPreproc ();

    if (!pCPreproc)
        return kMemError;

    return kNoError;
}

Error_t CPreproc::destroy( CPreproc*& pCPreproc )
{
    if (!pCPreproc)
        return kNoError;

    delete pCPreproc;
    pCPreproc   = 0;

    return kNoError;
}

Error_t CPreproc::init (const CAudioInfo* phAudioInfo)
{
    Error_t  rErr = kNoError;

    // sanity check
    if (phAudioInfo == 0)
        return kFunctionInvalidArgsError;

    // clean up
    reset();

    m_phAudioInfo       = phAudioInfo;
    m_bIsInitialized    = true;

    return rErr;
}

Error_t CPreproc::reset()
{
    
    m_bIsInitialized    = false;

    for (int i = 0; i < kNumPreprocSteps; i++)
        m_abIsStepActive[i] = false;

    m_phAudioInfo       = 0;

    return kNoError;

}

Error_t CPreproc::process( float **ppfInputBuffer, float **ppfOutputBuffer, int iNumOfFrames )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (!ppfInputBuffer || !ppfInputBuffer[0] || !ppfOutputBuffer || !ppfOutputBuffer[0] || iNumOfFrames < 0)
        return kFunctionInvalidArgsError;

    if (m_abIsStepActive[kPpDownmix])
    {
        downmix(ppfInputBuffer, ppfOutputBuffer[0], iNumOfFrames);
    }

    if (m_abIsStepActive[kPpNormalize])
    {
        if (m_abIsStepActive[kPpDownmix])
            normalize(ppfOutputBuffer[0], ppfOutputBuffer[0], iNumOfFrames);
        else
        {
            int iNumChannels = m_phAudioInfo->getNumChannels();
            for (int c = 0; c < iNumChannels; c++)
                normalize(ppfInputBuffer[c], ppfOutputBuffer[c], iNumOfFrames);
        }
    }

    return kNoError;
}

Error_t CPreproc::setStepActive( PreprocSteps_t eStep, bool bIsEnabled /*= true*/ )
{
    m_abIsStepActive[eStep] = bIsEnabled;

    return kNoError;
}

bool CPreproc::getStepActive( PreprocSteps_t eStep ) const
{
    return m_abIsStepActive[eStep];
}

Error_t CPreproc::downmix( float **ppfInBuff, float *pfOutBuff, int iNumOfFrames )
{
    int iNumChannels = m_phAudioInfo->getNumChannels();

    assert(iNumChannels > 0);
    assert(iNumOfFrames > 0);
    assert(ppfInBuff != 0);
    assert(ppfInBuff[0] != 0);
    assert(pfOutBuff != 0);

    CUtil::copyBuff(pfOutBuff, ppfInBuff[0], iNumOfFrames);

    for (int c = 1; c < iNumChannels; c++)
        CUtil::addBuff (pfOutBuff, ppfInBuff[c], iNumOfFrames);

    CUtil::mulBuffC(pfOutBuff, 1.F/iNumChannels, iNumOfFrames);

    return kNoError;
}

Error_t CPreproc::normalize( float *pfInBuff, float *pfOutBuff, int iNumOfFrames )
{
    int     iNumChannels = m_phAudioInfo->getNumChannels();
    double  dMax = 0;

    assert(iNumChannels > 0);
    assert(iNumOfFrames > 0);
    assert(pfInBuff != 0);
    assert(pfOutBuff != 0);

    if (m_abIsStepActive[kPpDownmix])
    {
        dMax = CUtil::getMax(pfInBuff,iNumOfFrames, true);
    }
    else
    {
        for (int c = 0; c < iNumChannels; c++)
        {
            double dTmp = 0;
            m_phAudioInfo->getResult(dTmp,CAudioInfo::kAbsMax, c);
            if (dMax < dTmp)
                dMax = dTmp;
        }
    }

    if (dMax > 0)
    {
        CUtil::copyBuff(pfOutBuff, pfInBuff, iNumOfFrames);
        CUtil::mulBuffC(pfOutBuff, static_cast<float>(1./dMax), iNumOfFrames);
    }

    return kNoError;
}

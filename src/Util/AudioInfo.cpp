#include <algorithm>
#include <cstdlib>

#include "Util.h"
#include "AudioInfo.h"

const std::string CAudioInfo::m_sResultNames[kNumInfoTypes] = {"Max", "Min", "AbsMax","Range","Mean","Rms","Std","LengthInS"};

std::string CAudioInfo::getResultName( InfoType_t eInfoType )
{
    return m_sResultNames[eInfoType];
}

CAudioInfo::CAudioInfo() :
    m_bIsInitialized(false),
    m_ppadResult(0),
    m_iNumChannels(0),
    m_fSampleRate(0)
{
    resetInstance ();
}

Error_t CAudioInfo::createInstance( CAudioInfo*& pCAudioInfo )
{
    pCAudioInfo = new CAudioInfo ();

    if (!pCAudioInfo)
        return kMemError;

    return kNoError;
}

Error_t CAudioInfo::destroyInstance( CAudioInfo*& pCAudioInfo )
{
    if (!pCAudioInfo)
        return kNoError;

    delete pCAudioInfo;
    pCAudioInfo   = 0;

    return kNoError;
}

Error_t CAudioInfo::initInstance (float fSampleRate, int iNumChannels)
{
    Error_t  rErr = kNoError;

    // sanity check
    if (fSampleRate <= 0 || iNumChannels <= 0)
        return kFunctionInvalidArgsError;

    // clean up
    resetInstance();

    m_fSampleRate   = fSampleRate;
    m_iNumChannels  = iNumChannels;
    m_ppadResult    = new double* [iNumChannels];
    for (int c = 0; c < iNumChannels; c++)
    {
        m_ppadResult[c] = new double [kNumInfoTypes];
        CUtil::setZero(m_ppadResult[c], kNumInfoTypes);
    }

    m_bIsInitialized    = true;

    return rErr;
}

Error_t CAudioInfo::resetInstance()
{
    
    m_bIsInitialized    = false;

    for (int c = 0; c < m_iNumChannels; c++)
    {
        delete [] m_ppadResult[c];
        m_ppadResult[c] = 0;
    }
    delete [] m_ppadResult;
    m_ppadResult    = 0;
    m_iNumChannels  = 0;

    return kNoError;

}

Error_t CAudioInfo::process( float **ppfAudio, long long int iLength )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (!ppfAudio || !ppfAudio[0] || iLength < 0)
        return kFunctionInvalidArgsError;

    for (int c = 0; c < m_iNumChannels; c++)
    {
        m_ppadResult[c][kMax]       = CUtil::getMax(ppfAudio[c], iLength);
        m_ppadResult[c][kMin]       = CUtil::getMin(ppfAudio[c], iLength);
        m_ppadResult[c][kAbsMax]    = std::max(std::abs(m_ppadResult[c][kMax]), std::abs(m_ppadResult[c][kMin]));
        m_ppadResult[c][kRange]     = std::abs(m_ppadResult[c][kMax]-m_ppadResult[c][kMin]);
        m_ppadResult[c][kMean]      = CUtil::getMean(ppfAudio[c], iLength);
        m_ppadResult[c][kStd]       = CUtil::getStd(ppfAudio[c], iLength);
        m_ppadResult[c][kRms]       = CUtil::getRms(ppfAudio[c], iLength);
        m_ppadResult[c][kLengthInS] = iLength/m_fSampleRate;
    }

    return kNoError;
}

Error_t CAudioInfo::getResult( double &dResult, InfoType_t eInfoType, int iChannelIdx /*= 0*/ )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (iChannelIdx < 0 || iChannelIdx > m_iNumChannels)
        return kFunctionInvalidArgsError;

    dResult = m_ppadResult[iChannelIdx][eInfoType];

    return kNoError;
}

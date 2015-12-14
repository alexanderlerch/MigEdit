#include <algorithm>


#include "Util.h"
#include "Fft.h"

#include "Ccf.h"



CCcf::CCcf () :
    m_pCFft (0),
    m_bNormalize(false),
    m_iFftLength(0)
{
    for (int i = 0; i < 2; i++)
    {
        m_aiLength[i]   = 0;
        m_apfData[i]    = 0;
    }
    reset();
}


CCcf::~CCcf ()
{
    reset();
}

Error_t CCcf::create (CCcf*& pCCcf)
{
    pCCcf   = 0;
    pCCcf   = new CCcf ();

    if (!pCCcf)
        return kMemError;

    return kNoError;
}


Error_t CCcf::destroy (CCcf*& pCCcf)
{
    delete pCCcf;
    pCCcf = 0;

    return kNoError;
}

Error_t CCcf::init( const int aiLength[2], bool bNormalize /*= false*/ )
{
    if (aiLength[0] <= 0 || aiLength[1] <= 0)
        return kFunctionInvalidArgsError;

    reset();

    m_iFftLength    = CUtil::nextPowOf2 ((std::max(aiLength[0], aiLength[1]))<<1);

    CFft::create(m_pCFft);
    m_pCFft->init(m_iFftLength, 1, CFft::kWindowHann, CFft::kNoWindow);

    for (int i = 0; i < 2; i++)
        m_apfData[i]    = new float [m_iFftLength];

    CUtil::copyBuff(m_aiLength, aiLength, 2);

    m_bNormalize    = bNormalize;

    return kNoError;
}

Error_t CCcf::reset()
{
    CFft::destroy(m_pCFft);

    for (int i = 0; i < 2; i++)
    {
        delete [] m_apfData[i];
        m_apfData[i]    = 0;
        m_aiLength[i]   = 0;
    }

    m_iFftLength    = 0;
    m_bNormalize    = false;

    return kNoError;
}

Error_t CCcf::process( const float *pfSrc1, const float *pfSrc2 )
{
    double afStd[2] = {0,0};
    assert ((pfSrc1 != 0) || (pfSrc2 != 0));

    // init internal  buffers
    for (int i = 0; i < 2; i++)
        CUtil::setZero(m_apfData[i], m_iFftLength);

    // copy data
    if (m_aiLength[0] >= m_aiLength[1])
    {
        CUtil::copyBuff(m_apfData[0], pfSrc1, m_aiLength[0]);
        CUtil::copyBuff(m_apfData[1], pfSrc2, m_aiLength[1]);
    }
    else
    {
        CUtil::copyBuff(m_apfData[0], pfSrc2, m_aiLength[1]);
        CUtil::copyBuff(m_apfData[1], pfSrc1, m_aiLength[0]);
    }

    if (m_bNormalize)
    {
        for (int i = 0; i < 2; i++)
            afStd[i] = CUtil::getStd(m_apfData[i], m_aiLength[i]);
    }
    // compute FFTs
    for (int i = 0; i < 2; i++)
        m_pCFft->doFft(m_apfData[i],m_apfData[i]);

    // conj comp mult
    m_pCFft->conjCompSpectrum(m_apfData[1], m_apfData[1]);
    m_pCFft->mulCompSpectrum (m_apfData[0], m_apfData[1]);

    // norm
    CUtil::mulBuffC(m_apfData[0], static_cast<float>(m_iFftLength), m_iFftLength);
    if (m_bNormalize && afStd[0]* afStd[1] > 0)
        CUtil::mulBuffC (m_apfData[0], 1.F/static_cast<float>(std::min(m_aiLength[0], m_aiLength[1])*afStd[0]*afStd[1]), m_iFftLength);

    // IFFT
    m_pCFft->doInvFft(m_apfData[0], m_apfData[0]);

    return kNoError;
}

Error_t CCcf::getResult( float *pfResult ) const
{
    int iCcfLength = getResultLength (m_aiLength);

    // copy result
    if (pfResult)
    {
        // Ccf begin starts at Fft length - shorter version length
        int iStartIdx       = m_iFftLength-std::min(m_aiLength[0],m_aiLength[1])+1,
            iStartLength    = std::min(iCcfLength,m_iFftLength - iStartIdx);

        CUtil::copyBuff(&pfResult[0], &m_apfData[0][iStartIdx], iStartLength);
        CUtil::copyBuff(&pfResult[iStartLength], m_apfData[0], (iCcfLength-iStartLength));

        if (m_aiLength[0] < m_aiLength[1])
            CUtil::flipBuff(pfResult, iCcfLength);
    }
    return kNoError;
}

int CCcf::getResultLength( const int aiLength[2] )
{
    return aiLength[0] + aiLength[1] - 1;
}

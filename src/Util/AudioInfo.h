#if !defined(__AudioInfo_hdr__)
#define __AudioInfo_hdr__

#include "ErrorDef.h"

class CAudioInfo
{
public:

    enum InfoType_t
    {
        kMax,                   //!< maximum sample value
        kMin,                   //!< minimum sample value
        kAbsMax,                //!< absolute maximum sample value
        kRange,                 //!< maximum range
        kMean,                  //!< average sample value
        kRms,                   //!< root mean square
        kStd,                   //!< standard deviation

        kNumInfoTypes
    };

    /*! creates a new AudioInfo instance
    \param CAudioInfo * & pCAudioInfo: pointer to the new instance
    \return Error_t
    */
    static Error_t createInstance (CAudioInfo*& pCAudioInfo);
    
    /*! destroys an AudioInfo instance
    \param CAudioInfo * & pCAudioInfo: pointer to the instance to be destroyed
    \return Error_t
    */
    static Error_t destroyInstance (CAudioInfo*& pCAudioInfo);
    
    /*! initializes an AudioInfo instance
    \param float fSamplerate: audio sample rate
    \param int iNumChannels: number of channels
    \return Error_t
    */
    Error_t initInstance (float fSampleRate, int iNumChannels);
    
    /*! resets an AudioInfo instance
    \return Error_t
    */
    Error_t resetInstance ();
 
    Error_t process (float **ppfAudio, int iLength);
    Error_t  getResult (double &dResult, InfoType_t eInfoType, int iChannelIdx = 0);

protected:
    CAudioInfo ();
    virtual ~CAudioInfo () {};

private:

    double  **m_ppadResult;

    int     m_iNumChannels;
    bool    m_bIsInitialized;
};

#endif // #if !defined(__AudioInfo_hdr__)




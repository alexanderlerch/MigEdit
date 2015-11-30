#if !defined(__AudioInfo_hdr__)
#define __AudioInfo_hdr__

#include <string>
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

        kLengthInS,             //!< length in samples

        kNumInfoTypes
    };

    /*! creates a new AudioInfo instance
    \param CAudioInfo * & pCAudioInfo: pointer to the new instance
    \return Error_t
    */
    static Error_t create (CAudioInfo*& pCAudioInfo);
    
    /*! destroys an AudioInfo instance
    \param CAudioInfo * & pCAudioInfo: pointer to the instance to be destroyed
    \return Error_t
    */
    static Error_t destroy (CAudioInfo*& pCAudioInfo);
    
    /*! initializes an AudioInfo instance
    \param float fSamplerate: audio sample rate
    \param int iNumChannels: number of channels
    \return Error_t
    */
    Error_t init (float fSampleRate, int iNumChannels);
    
    /*! resets an AudioInfo instance
    \return Error_t
    */
    Error_t reset ();
 
    Error_t process (float **ppfAudio, long long int iLength);

    int     getNumChannels() const;
    float   getSampleRate() const;
    Error_t getResult (double &dResult, InfoType_t eInfoType, int iChannelIdx = 0) const;
    std::string getResultName( InfoType_t eInfoType ) const;

protected:
    CAudioInfo ();
    virtual ~CAudioInfo () {};

private:
    CAudioInfo(const CAudioInfo& that);

    double  **m_ppadResult;

    int     m_iNumChannels;
    float   m_fSampleRate;
    bool    m_bIsInitialized;

    static const std::string m_sResultNames[kNumInfoTypes];
};

#endif // #if !defined(__AudioInfo_hdr__)




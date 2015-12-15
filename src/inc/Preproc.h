#if !defined(__Preproc_hdr__)
#define __Preproc_hdr__

#include "ErrorDef.h"

class CAudioInfo;

class CPreproc
{
public:
    /*! version number */
    enum PreprocSteps_t
    {
        kPpDownmix,                         
        kPpNormalize,                       

        kNumPreprocSteps
    };

    /*! creates a new CPreproc instance
    \param CPreproc * & pCInstance: pointer to the new instance
    \return Error_t
    */
    static Error_t create (CPreproc*& pCInstance);
    /*! destroys a CPreproc instance
    \param CPreproc * & pCInstance: pointer to instance to be destroyed
    \return Error_t
    */
    static Error_t destroy (CPreproc*& pCInstance);
    
    /*! initializes the class
    \param const CAudioInfo * pAudioInfo: class handle for info on the audio signal
    \return Error_t
    */
    Error_t init (const CAudioInfo* pAudioInfo);
    /*! resets the instance
    \return Error_t
    */
    Error_t reset ();

    /*! enables or disables a processing step
    \param PreprocSteps_t eStep: step to be set
    \param bool bIsEnabled
    \return Error_t
    */
    Error_t setStepActive (PreprocSteps_t eStep, bool bIsEnabled = true);
    /*! returns the state of processing steps
    \param PreprocSteps_t eStep: step to be returned the info for
    \return bool
    */
    bool    getStepActive (PreprocSteps_t eStep) const;

    /*! process the audio
    \param float * * ppfInputBuffer
    \param float * * ppfOutputBuffer
    \param int iNumOfFrames
    \return Error_t
    */
    Error_t process (float **ppfInputBuffer, float **ppfOutputBuffer, int iNumOfFrames);

protected:
    CPreproc ();
    virtual ~CPreproc (){};
private:
    Error_t normalize (float *pfInBuff, float *pfOutBuff, int iNumOfFrames);
    Error_t downmix (float **ppfInBuff, float *pfOutBuff, int iNumOfFrames);

    const CAudioInfo *m_phAudioInfo;
    bool    m_bIsInitialized;

    bool    m_abIsStepActive[kNumPreprocSteps];

};

#endif // #if !defined(__Preproc_hdr__)




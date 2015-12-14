#if !defined(__Ccf_hdr__)
#define __Ccf_hdr__

#include "ErrorDef.h"

class CFft;

class CCcf
{
public:
    static Error_t create (CCcf*& pCCcf);
    static Error_t destroy (CCcf*& pCCcf);
    Error_t init(const int aiLength[2], bool bNormalize = false);
    Error_t reset();
    Error_t process (const float *pfSrc1, const float *pfSrc2);

    static int  getResultLength (const int aiLength[2]);
    Error_t getResult (float *pfResult) const;
    

protected:
    CCcf ();
    virtual ~CCcf ();
private:
    int     m_aiLength[2];
    int     m_iFftLength;
    CFft    *m_pCFft;
    float   *m_apfData[2];

    bool    m_bNormalize;
};


#endif // #if !defined(__Ccf_hdr__)




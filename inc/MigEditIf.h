#if !defined(__MigEditIf_hdr__)
#define __MigEditIf_hdr__

#include "ErrorDef.h"

class CMigEditIf
{
public:
    /*! version number */
    enum Version_t
    {
        kMajor,                         //!< major version number
        kMinor,                         //!< minor version number
        kPatch,                         //!< patch version number

        kNumVersionInts
    };

    static const int  getVersion (const Version_t eVersionIdx);
    static const char* getBuildDate ();

    static Error_t createInstance (CMigEditIf*& pCInstance);
    static Error_t destroyInstance (CMigEditIf*& pCInstance);
    
    virtual Error_t initInstance (/*enter parameters here*/) = 0;
    virtual Error_t resetInstance () = 0;
    
    //virtual Error_t process (float **ppfInputBuffer, float **ppfOutputBuffer, int iNumberOfFrames) = 0;

protected:
    CMigEditIf () {};
    virtual ~CMigEditIf () {};
};

#endif // #if !defined(__MigEditIf_hdr__)



